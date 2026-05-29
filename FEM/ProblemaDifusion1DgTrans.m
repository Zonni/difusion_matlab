% CLASE: ProblemaDifusion1DGTrans
% Resuelve el problema de difusion neutronica 1D transitorio con G grupos de energia
% heredando de ProblemaDifusion1DG
%
% La clase proporciona:
%   - Ensamblaje de matrices de perdidas (L), fision con chi (M_chi)
%     y fision sin chi (M_fis)
%   - Rutina para integrador ODE (rhs_edo, solver) para ode15s
%   - Metodos para el esquema semi-implicito (4.14) (esqPC_2 para 2 grupos) 
%     que trabaja con PC = P*C
%   - Metodo para el esquema semi-implicito (53)-(54) (esqPC para G grupos) 
%     que trabaja con PC = P*C
%
% Notacion:
%   P         = matriz de masa del FEM
%   PC        = P * C  para los precursores
%   M_chi     = matriz de fision con chi.
%   M_fis{g}  = matriz de fision del grupo g sin chi
%   L         = matriz de perdidas (difusion + absorcion + out-scattering)
%               mas ganancias por in-scattering (bloques fuera de la diagonal)
%
% Se asume:
%   - Los neutrones inician en los grupos indicados por chi
%   - Scattering definido por la matriz sigma_s(g, h, mat) (de h a g)

classdef ProblemaDifusion1DgTrans < ProblemaDifusion1Dg
    properties
        nombre_banco        % Nombre del banco de datos
        masas_mat           % Matriz de masas P
        v                   % Vector de velocidades (G x 1)
        beta                % Vector de fracciones de neutrones diferidos
        lambda              % Vector de constantes de decaimiento
        n_precursor         % Numero de grupos de precursores
        chi                 % Chi (G x 1)

        % Funciones de seccion eficaz dependientes del tiempo (cell arrays {g, mat})
        D_func
        sigma_a_func
        nu_sigma_f_func
        sigma_s_func        % cell G x G x n_mat (scattering de h a g almacenado en {g, h, mat})

        % Resultados almacenados
        keff0
        Phi_t
        Pot_num
    end

    methods
        function obj = ProblemaDifusion1DgTrans(malla, materiales, elemento)
            % Hereda malla, materiales y elemento de la clase ProblemaDifusion1DG
            obj@ProblemaDifusion1Dg(malla, materiales, elemento);
            obj.chi = materiales.chi;
            obj.D_func = {};
            obj.sigma_a_func = {};
            obj.nu_sigma_f_func = {};
            obj.sigma_s_func = {};
        end

        function obj = ensamblar_matrices_cte(obj)
            % Ensambla las matrices que NO dependen del tiempo:
            % masas_mat = P, la matriz de masa del FEM
            % Es comun para todos los grupos y para los precursores
            %
            % La matriz P aparece en los terminos (1 / v_g) * P y en
            % los productos P*C_k de las ecuaciones de precursores
            %
            % Se aplican condiciones de contorno Dirichlet

            N = obj.malla.N;
            grado_l = obj.malla.grado_l;
            tam_celdas = obj.malla.tamano_celdas;
            Jac = tam_celdas / 2;

            Lag_eval = double(subs(obj.elemento.Lag, sym('x'), obj.elemento.xx'));

            M = N * grado_l + 1;

            masas = sparse(M, M);
            for e = 1 : N
                idx = obj.malla.nodos(e, :);
                masas_l = zeros(grado_l + 1);
                for i = 1 : grado_l + 1
                    for j = 1 : grado_l + 1
                        for k = 1 : grado_l + 1
                            w = obj.elemento.ww(k);
                            Li = Lag_eval(i, k);
                            Lj = Lag_eval(j, k);
                            masas_l(i,j) = masas_l(i,j) + w * Li * Lj;
                        end
                    end
                end
                masas_l = masas_l * Jac(e);
                masas(idx, idx) = masas(idx, idx) + masas_l;
            end

            % Aplicar condiciones de contorno
            nodos_borde = [1, M];
            for i = nodos_borde
                masas(i, :) = 0;
                masas(:, i) = 0;
                masas(i, i) = 1;
            end

            obj.masas_mat = masas;
        end

        function [L, M, M_fis] = ensamblar_matrices_t(obj, t)
            % Ensambla las matrices dependientes del tiempo en el instante t
            %
            % Salidas:
            %   L     : matriz de perdidas (G*M_g) x (G*M_g)
            %           Contiene los terminos de difusion (Ke), absorcion (Me) y
            %           out-scattering (So) en la diagonal de cada grupo.
            %           Fuera de la diagonal tiene los terminos de in-scattering ( - Ts )
            %           Corresponde a T^{h+1}
            %
            %   M     : matriz de fision con chi (G*M_g) x (G*M_g)
            %           Cada bloque (g, h) se construye como chi(g) * (nuSigma_f)_h * P
            %           Corresponde a F^{h+1}
            %
            %   M_fis : cell array de G matrices (M_g x M_g) 
            %           Se construyen sin chi como (nuSigma_f)_g * P
            %           Se usan para calcular la tasa total de fisiones (termino) y la potencia
            %
            % La matriz L incluye condiciones de contorno Dirichlet

            N = obj.malla.N;
            grado_l = obj.malla.grado_l;
            tam_celdas = obj.malla.tamano_celdas;
            ind_mat = obj.materiales.ind_materiales;
            Jac = tam_celdas / 2;

            Lag_eval = double(subs(obj.elemento.Lag, sym('x'), obj.elemento.xx'));
            dLag_eval = double(subs(obj.elemento.dLag, sym('x'), obj.elemento.xx'));

            G = obj.G;
            M_g = size(obj.masas_mat, 1);

            L = sparse(G * M_g, G * M_g);
            M = sparse(G * M_g, G * M_g);
            M_fis = cell(1, G);
            for g = 1 : G
                M_fis{g} = sparse(M_g, M_g);
            end

            for e = 1 : N
                idx = obj.malla.nodos(e, :);
                mat = ind_mat(e);

                % Obtener secciones eficaces en el instante t
                D_t = zeros(G, 1);
                sa_t = zeros(G, 1);
                nf_t = zeros(G, 1);
                for g = 1 : G
                    D_t(g) = obj.D_func{g, mat}(t);
                    sa_t(g) = obj.sigma_a_func{g, mat}(t);
                    nf_t(g) = obj.nu_sigma_f_func{g, mat}(t);
                end

                % Scattering de h a g: ss_t(g,h) = Sigma_s^{h -> g}
                ss_t = zeros(G, G);
                for g = 1 : G
                    for h = 1 : G
                        if g ~= h && ~isempty(obj.sigma_s_func{g, h, mat})
                            ss_t(g,h) = obj.sigma_s_func{g, h, mat}(t);
                        end
                    end
                end

                % Matrices locales para cada grupo
                for g = 1 : G
                    Ke = zeros(grado_l + 1);   % difusion
                    Me = zeros(grado_l + 1);   % absorcion
                    So = zeros(grado_l + 1);   % out-scattering
                    for i = 1 : grado_l + 1
                        for j = 1 : grado_l + 1
                            for k = 1 : grado_l + 1
                                w = obj.elemento.ww(k);
                                Li = Lag_eval(i, k);
                                Lj = Lag_eval(j, k);
                                dLi = dLag_eval(i, k);
                                dLj = dLag_eval(j, k);

                                Ke(i,j) = Ke(i,j) + w * D_t(g) * dLi * dLj;
                                Me(i,j) = Me(i,j) + w * sa_t(g) * Li * Lj;
                            end
                        end
                    end

                    % Out-scattering del grupo g
                    % ss_t(h,g) es el scattering de g a h
                    for h = 1 : G
                        if h ~= g && ss_t(h,g) ~= 0   
                            for i = 1 : grado_l + 1
                                for j = 1 : grado_l + 1
                                    for k = 1 : grado_l + 1
                                        w = obj.elemento.ww(k);
                                        So(i,j) = So(i,j) + w * ss_t(h,g) * Lag_eval(i,k) * Lag_eval(j,k);
                                    end
                                end
                            end
                        end
                    end

                    Ke = Ke / Jac(e);
                    Me = Me * Jac(e);
                    So = So * Jac(e);

                    L((g-1)*M_g+idx, (g-1)*M_g+idx) = L((g-1)*M_g+idx, (g-1)*M_g+idx) + Ke + Me + So;
                end

                % In-scattering en el grupo g desde h
                % ss_t(g,h) es scattering de h a g
                for g = 1 : G
                    for h = 1 : G
                        if g ~= h && ss_t(g,h) ~= 0
                            Ts = zeros(grado_l + 1);
                            for i = 1 : grado_l + 1
                                for j = 1 : grado_l + 1
                                    for k = 1 : grado_l + 1
                                        w = obj.elemento.ww(k);
                                        Ts(i,j) = Ts(i,j) + w * ss_t(g,h) * Lag_eval(i,k) * Lag_eval(j,k);
                                    end
                                end
                            end
                            Ts = Ts * Jac(e);
                            L((g - 1) * M_g + idx, (h - 1) * M_g + idx) = L((g - 1) * M_g + idx, (h - 1) * M_g + idx) - Ts;
                        end
                    end
                end

                % Fision con chi: M_{g,h} = chi_g * (nuSigma_f)_h * P
                for g = 1 : G
                    chi_g = obj.chi(g);
                    if chi_g ~= 0
                        for h = 1 : G
                            nf_h = nf_t(h);
                            if nf_h ~= 0
                                Te = zeros(grado_l + 1);
                                for i = 1 : grado_l + 1
                                    for j = 1 : grado_l + 1
                                        for k = 1 : grado_l + 1
                                            w = obj.elemento.ww(k);
                                            Te(i,j) = Te(i,j) + w * chi_g * nf_h * Lag_eval(i,k) * Lag_eval(j,k);
                                        end
                                    end
                                end
                                Te = Te * Jac(e);
                                M((g-1)*M_g+idx, (h-1)*M_g+idx) = M((g-1)*M_g+idx, (h-1)*M_g+idx) + Te;
                            end
                        end
                    end
                end

                % Fision sin chi: (nuSigma_f)_g * P
                for g = 1 : G
                    nf_g = nf_t(g);
                    if nf_g ~= 0
                        Te = zeros(grado_l + 1);
                        for i = 1 : grado_l + 1
                            for j = 1 : grado_l + 1
                                for k = 1 : grado_l + 1
                                    w = obj.elemento.ww(k);
                                    Te(i,j) = Te(i,j) + w * nf_g * Lag_eval(i,k) * Lag_eval(j,k);
                                end
                            end
                        end
                        Te = Te * Jac(e);
                        M_fis{g}(idx, idx) = M_fis{g}(idx, idx) + Te;
                    end
                end
            end

            % Aplicar condiciones de contorno Dirichlet
            nodos_borde = [1, M_g];
            for g = 1 : G
                base = (g-1)*M_g;
                for i = nodos_borde
                    L(base+i, :) = 0;
                    L(base+i, base+i) = 1;
                    for h = 1 : G
                        col_inicio = (h-1)*M_g + 1;
                        col_fin    = col_inicio + M_g - 1;
                        M(base+i, col_inicio:col_fin) = 0;
                    end
                end
            end
            for g = 1 : G
                M_fis{g}(nodos_borde, :) = 0;
            end
        end

        function y0 = cond_ini(obj, Phi_est, t0)
            % Calcula la condicion inicial para el integrador ODE (rhs_edo)
            % Devuelve [Phi_est; C0(:)]

            G = obj.G;
            n_dof = size(obj.masas_mat, 1);
            K = obj.n_precursor;

            if K == 0
                y0 = Phi_est;
            else
                % Ensamblar matrices en t0 para obtener el termino de fision
                [~, ~, M_fis] = obj.ensamblar_matrices_t(t0);
                term = zeros(n_dof, 1);
                for g = 1 : G
                    term = term + M_fis{g} * Phi_est((g - 1) * n_dof + 1 : g * n_dof);
                end

                C0 = zeros(n_dof, K);
                for k = 1 : K
                    rhs = (obj.beta(k) / obj.lambda(k)) * term;
                    C0(:, k) = obj.masas_mat \ rhs;
                end
                y0 = [Phi_est; C0(:)];
            end
        end

        function dydt = rhs_edo(obj, t, y)
            % Calcula el lado derecho del sistema para el integrador temporal
            %   dPhi/dt = V^{-1} * (M_chi * Phi + sum_k lambda_k * chi * P * C_k - L * Phi)
            %   dC_k/dt = beta_k * P^{-1} * F_1 * Phi - lambda_k * C_k
            % donde V = diag(1 / v_g) * P y los precursores se incorporan
            % como chi(g) * (P * C * lambda) en cada grupo

            G = obj.G;
            n_dof = size(obj.masas_mat, 1);
            K = obj.n_precursor;

            Phi = y(1 : G * n_dof);
            if K > 0
                C = reshape(y(G * n_dof + 1 : end), n_dof, K);
            else
                C = [];
            end

            [L, ~, M_fis] = obj.ensamblar_matrices_t(t);

            % Tasa total de fisiones (sin chi)
            termino = zeros(n_dof, 1);
            for g = 1 : G
                termino = termino + M_fis{g} * Phi((g - 1) * n_dof + 1 : g * n_dof);
            end

            % Fuente de neutrones retardados: P * C * lambda
            if K > 0
                term_retardados = obj.masas_mat * (C * obj.lambda(:));
            else
                term_retardados = zeros(n_dof, 1);
            end

            beta_total = sum(obj.beta);

            dPhi = zeros(G * n_dof, 1);
            for g = 1 : G
                idx = (g - 1) * n_dof + 1 : g * n_dof;
                perdidas = L(idx, :) * Phi;
                %   (1 - beta_total) * chi(g) * termino + chi(g) * term_retardados
                term_g = (1 - beta_total) * obj.chi(g) * termino + obj.chi(g) * term_retardados;
                % dPhi/dt = v_g * P^{-1} * ( - perdidas + fuente )
                dPhi(idx) = obj.v(g) * (obj.masas_mat \ ( - perdidas + term_g ));
            end

            if K > 0
                cociente = obj.masas_mat \ termino;
                dCdt = zeros(n_dof, K);
                for k = 1 : K
                    dCdt(:, k) = obj.beta(k) * cociente - obj.lambda(k) * C(:, k);
                end
                dydt = [dPhi; dCdt(:)];
            else
                dydt = dPhi;
            end
        end

        function [t, y] = solver(obj, int_t, y0, opciones)
            if nargin < 4
                opciones = odeset();
            end
            [t, y] = ode15s(@(t,y) obj.rhs_edo(t, y), int_t, y0, opciones);
        end

        function P_norm = pot_norm(obj, t, y)
            % Calcula la potencia normalizada P_norm = P(t) / P(0)
            % Integra sum_g (nuSigma_f)_g * phi_g usando cuadratura de Gauss

            G = obj.G;
            n_dof = size(obj.masas_mat, 1);
            ind_mat = obj.materiales.ind_materiales;

            xx = obj.elemento.xx;
            ww = obj.elemento.ww;
            Lag_eval = double(subs(obj.elemento.Lag, sym('x'), xx'));

            P_t = zeros(length(t), 1);
            for i = 1 : length(t)
                integral = 0;
                for e = 1 : obj.malla.N
                    mat = ind_mat(e);
                    h_e = obj.malla.tamano_celdas(e);
                    Jac = h_e / 2;
                    f_sum = zeros(size(xx));
                    for g = 1 : G
                        nf_t = obj.nu_sigma_f_func{g, mat}(t(i));
                        phi_g = y(i, (g - 1) * n_dof + 1 : g * n_dof)';
                        idx = obj.malla.nodos(e, :);
                        f_g = Lag_eval' * phi_g(idx);
                        f_sum = f_sum + nf_t * f_g;
                    end
                    integral = integral + sum(ww .* f_sum) * Jac;
                end
                P_t(i) = integral;
            end
            P_norm = P_t / P_t(1);
        end

        function graficar_t(obj, t, y)
            % Grafica potencia normalizada y flujo en varios instantes
            % Opcional: guarda la figura en formato EPS

            G = obj.G;
            n_dof = size(obj.masas_mat, 1);
            P_norm = obj.pot_norm(t, y);

            % Instantes de tiempo
            tiempos_plot = linspace(t(1), t(end), 5);
            idx_plot = zeros(size(tiempos_plot));
            for i = 1 : length(tiempos_plot)
                [~, idx_plot(i)] = min(abs(t - tiempos_plot(i)));
            end

            x_nodos = obj.malla.x_nodos;

            figure('Name', ['Resultados transitorio 1D', num2str(G), 'G: ', obj.nombre_banco]);

            subplot(1, 2, 1);
            plot(t, P_norm, 'b-', 'LineWidth', 2);
            xlabel('Tiempo (s)');
            ylabel('Potencia normalizada');
            grid on;

            subplot(1, 2, 2);
            hold on;
            colores = jet(length(idx_plot));
            grupos_mostrar = unique([1, G]);
            estilos = {'-', '--', ':'};
            for ig = 1 : length(grupos_mostrar)
                g = grupos_mostrar(ig);
                for i = 1 : length(idx_plot)
                    idx = idx_plot(i);
                    flujo_g = y(idx, (g - 1) * n_dof + 1 : g * n_dof);
                    plot(x_nodos, flujo_g, 'Color', colores(i,:), ...
                         'LineStyle', estilos{min(ig, length(estilos))}, ...
                         'DisplayName', sprintf('G%d t=%.3f s', g, t(idx)), 'LineWidth', 2);
                end
            end
            xlabel('x (cm)');
            ylabel('Flujo');
            legend('Location', 'best');
            grid on;

            saveas(gcf, fullfile('results/figures', ['Resultados_', obj.nombre_banco, '.eps']), 'epsc'); % Opcional
        end
        
        % --- Metodos para el esquema semi‑implicito (4.14) para 2 grupos ---

        function y0 = cond_ini_PC0_2(obj, phi1_0, phi2_0, t0)
            % Devuelve el estado inicial con PC = P*C para el esquema de 2 grupos
            % Se obtiene a partir de cond_ini (que devuelve C) y se multiplica por P
            n_dof = size(obj.masas_mat, 1);
            K = obj.n_precursor;
            if K == 0
                y0 = [phi1_0; phi2_0];
            else
                Phi0 = [phi1_0; phi2_0];
                y0_C = obj.cond_ini(Phi0, t0);
                C0_mat = reshape(y0_C(2 * n_dof + 1 : end), n_dof, K);
                PC0 = obj.masas_mat * C0_mat;   
                y0 = [phi1_0; phi2_0; PC0(:)];
            end
        end
        
        function [phi1_act, phi2_act, PC_act] = esqPC_2(obj, paso_t, phi1_n, phi2_n, PC_n, t_act)
            % Esquema semi‑implicito para 2 grupos (4.14)
            %
            % Utiliza las matrices L y M obtenidas de ensamblar_matrices_t y
            % construye T^{n+1} y E^n y se actualiza PC
            %
            % Entradas:
            %   paso_t : Delta t
            %   phi1_n, phi2_n : flujos en t^n
            %   PC_n   : P*C en t^n (n_dof x K)
            %   t_act  : t^{n+1}
            % Salidas:
            %   phi1_act, phi2_act : flujos en t^{n+1}
            %   PC_act : P*C en t^{n+1}

            n_dof = size(obj.masas_mat, 1);
            K = obj.n_precursor;
            P = obj.masas_mat;
        
            % Matrices en t_act
            [L, M, ~] = obj.ensamblar_matrices_t(t_act);
        
            % Extraer bloques 2x2
            L11 = L(1 : n_dof, 1 : n_dof);
            L21 = L(n_dof + 1 : 2 * n_dof, 1 : n_dof);
            L22 = L(n_dof + 1 : 2 * n_dof, n_dof + 1 : 2 * n_dof);
            M11 = M(1 : n_dof, 1 : n_dof);
            M12 = M(1 : n_dof, n_dof + 1 : 2 * n_dof);
        
            % Coeficiente a_hat de (4.14)
            if K > 0
                a = 1 - sum(obj.beta) + sum(obj.beta .* (1 - exp(-obj.lambda * paso_t)));
            else
                a = 1;
            end
        
            % Velocidades de los dos grupos
            v1 = obj.v(1);
            v2 = obj.v(2);
        
            % Construir T^{n+1} = (1 / dt) * V_inv + L - a * M
            inv_v1_dt = 1 / (v1 * paso_t);
            inv_v2_dt = 1 / (v2 * paso_t);
        
            T11 = inv_v1_dt * P + L11 - a * M11;
            T12 = - a * M12;
            T21 = L21;
            T22 = inv_v2_dt * P + L22;
            T = [T11, T12; T21, T22];
        
            % Construir E^n = (1 / dt) * V_inv * Phi_n + precursores
            E1 = inv_v1_dt * (P * phi1_n);
            E2 = inv_v2_dt * (P * phi2_n);
        
            if K > 0
                % sum_p lambda_p * exp(-lambda_p * dt) * (P*C_p^n)
                prec_term = PC_n * (obj.lambda(:) .* exp(-obj.lambda(:) * paso_t));
                E1 = E1 + prec_term;   % solo grupo rapido (chi = [1; 0])
            end
        
            E = [E1; E2];
        
            % Resolver flujos
            Phi_act = T \ E;
            phi1_act = Phi_act(1 : n_dof);
            phi2_act = Phi_act(n_dof + 1 : end);
        
            % Actualizar PC
            if K > 0
                termino = M11 * phi1_act + M12 * phi2_act;
                PC_act = zeros(n_dof, K);
                for p = 1 : K
                    PC_act(:, p) = exp(-obj.lambda(p) * paso_t) * PC_n(:, p) + ...
                                   (obj.beta(p) / obj.lambda(p)) * ...
                                   (1 - exp(-obj.lambda(p) * paso_t)) * termino;
                end
            else
                PC_act = [];
            end
        end

        % --- Metodos para el esquema semi-implicito (53)-(54) para G grupos ---

        function y0 = cond_ini_PC0(obj, Phi_est, t0)
            % Devuelve el estado inicial con PC = P*C para el esquema de G grupos
            % Se obtiene a partir de cond_ini (que devuelve C) y se multiplica por P
            G = obj.G;
            n_dof = size(obj.masas_mat, 1);
            K = obj.n_precursor;
            if K == 0
                y0 = Phi_est;
            else
                y0_C = obj.cond_ini(Phi_est, t0);
                C0_mat = reshape(y0_C(G * n_dof + 1 : end), n_dof, K);
                PC0 = obj.masas_mat * C0_mat;  
                y0 = [Phi_est; PC0(:)];
            end
        end

        function [Phi_act, PC_act] = esqPC(obj, paso_t, Phi_n, PC_n, t_act)
            % Esquema semi‑implicito de ecuaciones (53) y (54)
            %
            % Entradas:
            %   paso_t : Delta t
            %   Phi_n  : flujo en t^n (vector G*n_dof x 1)
            %   PC_n   : P*C en t^n (n_dof x K)
            %   t_act  : tiempo t^{n+1}
            %
            % Salidas:
            %   Phi_act : flujo en t^{n+1}
            %   PC_act  : P*C en t^{n+1} (n_dof x K)

            G = obj.G;
            n_dof = size(obj.masas_mat, 1);
            K = obj.n_precursor;
            P = obj.masas_mat;

            % Matrices en t^{n+1}
            [L, M_chi, M_fis] = obj.ensamblar_matrices_t(t_act);

            % Ecuacion (53): Sistema lineal para el flujo
            % T = (1 / dt) * V + L - M_chi
            V_inv_dt = sparse(G * n_dof, G * n_dof);
            for g = 1:G
                V_inv_dt((g - 1) * n_dof + 1 : g * n_dof, (g - 1) * n_dof + 1 : g * n_dof) = ...
                    (1 / (obj.v(g) * paso_t)) * P;
            end
            T = V_inv_dt + L - M_chi;

            % (1 / dt) * V * Phi_n + sum_k lambda_k * chi * (P * C_k^n)
            E = V_inv_dt * Phi_n;
            if K > 0
                contrib_prec = PC_n * obj.lambda(:);   % n_dof x 1
                for g = 1 : G
                    if obj.chi(g) ~= 0
                        idx_g = (g - 1) * n_dof + 1 : g * n_dof;
                        E(idx_g) = E(idx_g) + obj.chi(g) * contrib_prec;
                    end
                end
            end
            
            % Resolver flujo
            Phi_act = T \ E;

            % Ecuacion (54): Actualizacion de PC
            % ( (1 / dt) + lambda_k ) * (P*C_k^{h+1}) = (1 / dt) * (P*C_k^h) + beta_k * termino
            if K > 0
                termino = zeros(n_dof, 1);
                for g = 1 : G
                    idx_g = (g - 1) * n_dof + 1 : g * n_dof;
                    termino = termino + M_fis{g} * Phi_act(idx_g);
                end

                PC_act = zeros(n_dof, K);
                for k = 1 : K
                    alpha = 1 / paso_t + obj.lambda(k);
                    PC_act(:, k) = (1 / alpha) * ( (1 / paso_t) * PC_n(:, k) + obj.beta(k) * termino );
                end
            else
                PC_act = [];
            end
        end
    end
end