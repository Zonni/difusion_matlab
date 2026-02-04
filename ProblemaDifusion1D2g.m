% CLASE: ProblemaDifusion1D2g
% Resuelve la difusión de neutrones en 1D para 2 grupos de energía mediante dos métodos: Elementos Finitos (FEM) o 
% Diferencias Finitas (FDM).

classdef ProblemaDifusion1D2g
    properties
        malla
        materiales
        elemento
        A           
        B           
        phi_inc     
        keff        
    end
    
    methods
        % CONSTRUCTOR: Asigna malla, materiales y definición del elemento.
        function s = ProblemaDifusion1D2g(malla, materiales, elemento)
            s.malla = malla;
            s.materiales = materiales;
            s.elemento = elemento;
        end
        
        % MÉTODO FEM: Ensamble mediante Elementos Finitos.
        function s = ensamblar_matrices_fem(s)
            % 1. Extracción de parámetros de malla y materiales (vectores de 2 filas)
            N = s.malla.N;
            grado_l = s.malla.grado_l;
            tamano_celdas = s.malla.tamano_celdas;
            D1 = s.materiales.D(1,:);
            D2 = s.materiales.D(2,:);
            sigma_a1 = s.materiales.sigma_a(1,:);
            sigma_a2 = s.materiales.sigma_a(2,:);
            nu_sigma_f1 = s.materiales.nu_sigma_f(1,:);
            nu_sigma_f2 = s.materiales.nu_sigma_f(2,:);
            sigma_s12 = s.materiales.sigma_s12; 
            ind_materiales = s.materiales.ind_materiales;
            
            % 2. Evaluación de funciones de forma y Jacobiano
            Lag_eval = double(subs(s.elemento.Lag, sym('x'), s.elemento.xx'));
            dLag_eval = double(subs(s.elemento.dLag, sym('x'), s.elemento.xx'));
            Jac = tamano_celdas / 2;
            
            % 3. Prealocación de matrices
            M = N*grado_l + 1; % número de nodos por grupo
            A = sparse(2*M, 2*M);
            B = sparse(2*M, 2*M);
            
            % 4. Bucle de integración por elementos
            for e = 1:N
                idx = s.malla.nodos(e,:);
                % Inicialización de matrices locales
                Ke1 = zeros(grado_l+1); Ke2 = zeros(grado_l+1);
                Me1 = zeros(grado_l+1); Me2 = zeros(grado_l+1);
                Te1 = zeros(grado_l+1); Te2 = zeros(grado_l+1);
                Ts12 = zeros(grado_l+1); 
                
                % Cuadratura numérica
                for i = 1:grado_l+1
                    for j = 1:grado_l+1
                        for k = 1:grado_l+1
                            % Términos Grupo 1 
                            Ke1(i,j) = Ke1(i,j) + s.elemento.ww(k)*D1(ind_materiales(e))*dLag_eval(i,k)*dLag_eval(j,k);
                            Me1(i,j) = Me1(i,j) + s.elemento.ww(k)*(sigma_a1(ind_materiales(e))+sigma_s12(ind_materiales(e)))*Lag_eval(i,k)*Lag_eval(j,k);
                            Te1(i,j) = Te1(i,j) + s.elemento.ww(k)*nu_sigma_f1(ind_materiales(e))*Lag_eval(i,k)*Lag_eval(j,k);
                            
                            % Términos Grupo 2 
                            Ke2(i,j) = Ke2(i,j) + s.elemento.ww(k)*D2(ind_materiales(e))*dLag_eval(i,k)*dLag_eval(j,k);
                            Me2(i,j) = Me2(i,j) + s.elemento.ww(k)*sigma_a2(ind_materiales(e))*Lag_eval(i,k)*Lag_eval(j,k);
                            Te2(i,j) = Te2(i,j) + s.elemento.ww(k)*nu_sigma_f2(ind_materiales(e))*Lag_eval(i,k)*Lag_eval(j,k);
                            
                            % Término de Acoplamiento (Fuente G1 -> G2)
                            Ts12(i,j) = Ts12(i,j) + s.elemento.ww(k)*sigma_s12(ind_materiales(e))*Lag_eval(i,k)*Lag_eval(j,k);
                        end
                    end
                end
                
                % Normalización con Jacobiano
                Ke1 = Ke1 / Jac(e); Me1 = Me1 * Jac(e); Te1 = Te1 * Jac(e);
                Ke2 = Ke2 / Jac(e); Me2 = Me2 * Jac(e); Te2 = Te2 * Jac(e);
                Ts12 = Ts12 * Jac(e);
                
                % 5. Ensamblaje en matrices globales por bloques
                A(idx, idx) = A(idx, idx) + Ke1 + Me1;          
                A(M+idx, M+idx) = A(M+idx, M+idx) + Ke2 + Me2;   
                B(idx, idx) = B(idx, idx) + Te1;                 
                B(idx, M+idx) = B(idx, M+idx) + Te2;             
                A(M+idx, idx) = A(M+idx, idx) - Ts12;            
            end
            s.A = A;
            s.B = B;
        end
        
        % MÉTODO FDM: Ensamble mediante Diferencias Finitas.
        function s = ensamblar_matrices_fdm(s)
            % 1. Definición de variables espaciales y materiales
            N = s.malla.N;
            L = s.malla.L;
            delta_x = L / (N+2); 
            x_vec = linspace(0, L, N+2);
            tamano_celdas = s.malla.tamano_celdas;
            
            % Extracción de vectores de material
            D1 = s.materiales.D(1,:); D2 = s.materiales.D(2,:);
            sigma_a1 = s.materiales.sigma_a(1,:); sigma_a2 = s.materiales.sigma_a(2,:);
            nu_sigma_f1 = s.materiales.nu_sigma_f(1,:); nu_sigma_f2 = s.materiales.nu_sigma_f(2,:);
            sigma_s12 = s.materiales.sigma_s12;
            ind_materiales = s.materiales.ind_materiales;
            
            n_groups = 2;
            n_dofs = n_groups * N;
            
            % 2. Matrices dispersas
          
            A = spalloc(n_dofs, n_dofs, 3*n_dofs);
            B = spalloc(n_dofs, n_dofs, 3*n_dofs);
            
            % 3. Ecuaciones para Nodo 1
            % Grupo 1
            A(1, 1) = 2*D1(ind_materiales(1))/delta_x^2 + sigma_a1(ind_materiales(1)) + sigma_s12(ind_materiales(1));
            A(1, 2) = - D1(ind_materiales(1))/delta_x^2;
            B(1, 1) = nu_sigma_f1(ind_materiales(1));
            % Grupo 2 
            A(1+N, 1) = - sigma_s12(ind_materiales(1));
            A(1+N, 1+N) = 2*D2(ind_materiales(1))/delta_x^2 + sigma_a2(ind_materiales(1));
            A(1+N, 2+N) = - D2(ind_materiales(1))/delta_x^2;
            B(1+N, 1+N) = nu_sigma_f2(ind_materiales(1));
            
            % 4. Ecuaciones para Nodo N 
            % Grupo 1
            A(N, N) = 2*D1(ind_materiales(N))/delta_x^2 + sigma_a1(ind_materiales(N)) + sigma_s12(ind_materiales(N));
            A(N, N-1) = - D1(ind_materiales(N))/delta_x^2;
            B(N, N) = nu_sigma_f1(ind_materiales(N));
            % Grupo 2
            A(2*N, N) = - sigma_s12(ind_materiales(N));
            A(2*N, 2*N) = 2*D2(ind_materiales(N))/delta_x^2 + sigma_a2(ind_materiales(N));
            A(2*N, 2*N-1) = - D2(ind_materiales(N))/delta_x^2;
            B(2*N, 2*N) = nu_sigma_f2(ind_materiales(N));
            
            % 5. Bucle Principal (Nodos internos 2...N-1)
            for i = 2:N-1
                
                % Lógica manual para determinar material según posición x
                if (x_vec(i+1)< tamano_celdas(1)) || (x_vec(i+1)>sum(tamano_celdas(1:N-1)))
                    mat = 1;
                else
                    mat = 2;
                end
                
                % Determinar material del nodo anterior (mat_ant)
                if (x_vec(i+1)< tamano_celdas(1)) || (x_vec(i+1)>sum(tamano_celdas(1:N-1)))
                    mat_ant = 1;
                else
                    mat_ant = 2;
                end
                
                % Determinar material del nodo siguiente (mat_sig)
                if (x_vec(i+1)< tamano_celdas(1)) || (x_vec(i+1)>sum(tamano_celdas(1:N-1)))
                    mat_sig = 1;
                else
                    mat_sig = 2;
                end
                
                % Creacion de matrices usando promedios de difusión
                % A11: Ecuación Grupo 1
                A(i, i-1) = - 0.5 * (D1(ind_materiales(mat)) + D1(ind_materiales(mat_ant))) / delta_x^2 ;
                A(i, i) = 0.5 *(D1(ind_materiales(mat_ant)) + 2*D1(ind_materiales(mat)) + D1(ind_materiales(mat_sig)))/ delta_x^2 + sigma_a1(ind_materiales(mat)) + sigma_s12(ind_materiales(mat));
                A(i, i+1) = - 0.5 * (D1(ind_materiales(mat)) + D1(ind_materiales(mat_sig))) / delta_x^2 ;
                
                % A21: Término de Scattering G1 -> G2
                A(i+N, i) = - sigma_s12(ind_materiales(mat));
                
                % A22: Ecuación Grupo 2
                A(i+N, i-1+N) = - 0.5 * (D2(ind_materiales(mat)) + D2(ind_materiales(mat_ant))) / delta_x^2 ;
                A(i+N, i+N) = 0.5 *(D2(ind_materiales(mat_ant)) + 2*D2(ind_materiales(mat)) + D2_sig)/ delta_x^2 + sigma_a2(ind_materiales(mat));
                A(i+N, i+1+N) = - 0.5 * (D2(ind_materiales(mat)) + D2(ind_materiales(mat_sig))) / delta_x^2 ;
                
                % Matriz B: Fisión (Grupo 1 recibe aporte de ambos grupos)
                % B11
                B(i, i) = nu_sigma_f1(ind_materiales(mat));
                % B12
                B(i, i+N) = nu_sigma_f2(ind_materiales(mat));
                
                s.A = A;
                s.B = B;
            end
        end
        
        % MÉTODO CC: Aplica condiciones de contorno Dirichlet (Phi=0).
        function s = aplicar_cc(s)
            M = s.malla.N*s.malla.grado_l + 1;
            % Grupo 1
            s.A(1,:) = 0; s.A(1,1)=1; s.B(1,:)=0;
            s.A(M,:) = 0; s.A(M,M)=1; s.B(M,:)=0;
            
            % Grupo 2
            s.A(M+1,:) = 0; s.A(M+1,M+1)=1; s.B(M+1,:)=0;
            s.A(end,:) = 0; s.A(end,end)=1; s.B(end,:)=0;
        end
        
        % MÉTODO SOLVER: Resuelve autovalores y normaliza flujos.
        function s = resolver_autovalor(s,nm)
            [phi,K] = eigs(s.B, s.A, nm, 'largestreal');
            s.keff = diag(K);
            
            for i = 1:nm
                phi(:,i)=phi(:,i) / max(phi(:,i));
                if phi(2,i)<0
                    phi(:,i)=abs(phi(:,i));
                end
            end
            s.phi_inc = phi;
        end
        
        % MÉTODO GRAFICAR: Visualiza los flujos separados por grupos.
        function s = graficar(s,modos)
            M = s.malla.N*s.malla.grado_l + 1;
            x = s.malla.x_nodos;
            figure; hold on; grid on;
            colores_grupo = {'b','r'}; % b = grupo1, r = grupo2
            estilos_modo = {'-','--',':','-.'};
            
            for idx = 1:length(modos)
                modo = modos(idx);
                estilo = estilos_modo{mod(idx-1,length(estilos_modo))+1};
                % Gráfica Grupo 1
                plot(x, s.phi_inc(1:M,modo), 'Color', colores_grupo{1},'LineStyle', estilo, 'LineWidth', 1.5);
                % Gráfica Grupo 2
                plot(x, s.phi_inc(M+1:end,modo), 'Color', colores_grupo{2},'LineStyle', estilo, 'LineWidth', 1.5);
            end
            xlabel('x (cm)'); ylabel('\phi(x)');
            title('Modos FEM, ambos grupos');
            
            %leyenda dinámica
            leyenda = {};
            for idx = 1:length(modos)
                leyenda{end+1} = sprintf('Modo %d - Grupo 1', modos(idx));
                leyenda{end+1} = sprintf('Modo %d - Grupo 2', modos(idx));
            end
            legend(leyenda,'Location','best');
        end
    end
end