% Archivo: ProblemaDifusion1DG.m
% CLASE: ProblemaDifusion1DG
% Resuelve la difusion de neutrones en 1D para G grupos de energia mediante FEM

classdef ProblemaDifusion1Dg
    properties
        malla
        materiales
        elemento
        G
        A
        B
        phi_inc
        keff
    end

    methods
        function obj = ProblemaDifusion1Dg(malla, materiales, elemento)
            obj.malla = malla;
            obj.materiales = materiales;
            obj.elemento = elemento;
            obj.G = size(materiales.D, 1);
        end

        function obj = ensamblar_matrices_fem(obj)
            N = obj.malla.N;
            grado_l = obj.malla.grado_l;
            tam_celdas = obj.malla.tamano_celdas;
            G = obj.G;
            M_nodos = N * grado_l + 1;

            D = obj.materiales.D;
            sigma_a = obj.materiales.sigma_a;
            nu_sigma_f = obj.materiales.nu_sigma_f;
            sigma_s = obj.materiales.sigma_s;
            chi = obj.materiales.chi;
            ind_mat = obj.materiales.ind_materiales;

            Lag_eval = double(subs(obj.elemento.Lag, sym('x'), obj.elemento.xx'));
            dLag_eval = double(subs(obj.elemento.dLag, sym('x'), obj.elemento.xx'));
            Jac = tam_celdas / 2;

            A = sparse(G * M_nodos, G * M_nodos);
            B = sparse(G * M_nodos, G * M_nodos);

            for e = 1:N
                idx = obj.malla.nodos(e, :);
                mat = ind_mat(e);

                for g = 1:G
                    D_g = D(g, mat);
                    sa_g = sigma_a(g, mat);

                    Ke = zeros(grado_l+1);
                    Me = zeros(grado_l+1);
                    for i = 1:grado_l+1
                        for j = 1:grado_l+1
                            for k = 1:grado_l+1
                                w = obj.elemento.ww(k);
                                Ke(i,j) = Ke(i,j) + w * D_g * dLag_eval(i,k)*dLag_eval(j,k);
                                Me(i,j) = Me(i,j) + w * sa_g * Lag_eval(i,k)*Lag_eval(j,k);
                            end
                        end
                    end
                    Ke = Ke / Jac(e);
                    Me = Me * Jac(e);

                    % sigma_s(h,g,mat) es scattering de g a h (perdida para g)
                    scat_per = zeros(grado_l+1);
                    for h = 1 : G
                        if h ~= g
                            ss_gh = sigma_s(h, g, mat);  
                            if ss_gh ~= 0
                                for i = 1 : grado_l + 1
                                    for j = 1 : grado_l + 1
                                        for k = 1 : grado_l + 1
                                            w = obj.elemento.ww(k);
                                            scat_per(i,j) = scat_per(i,j) + w * ss_gh * Lag_eval(i,k) * Lag_eval(j,k);
                                        end
                                    end
                                end
                            end
                        end
                    end
                    scat_per = scat_per * Jac(e);
                    A_g_g = Ke + Me + scat_per;
                    A((g - 1) * M_nodos+idx, (g - 1) * M_nodos+idx) = ...
                        A((g - 1)*M_nodos+idx, (g - 1) * M_nodos+idx) + A_g_g;
                end

                % sigma_s(g,h,mat) es scattering de h a g (ganancia para g)
                for g = 1 : G
                    for h = 1 : G
                        if g ~= h
                            ss_hg = sigma_s(g, h, mat); 
                            if ss_hg ~= 0
                                Ts = zeros(grado_l + 1);
                                for i = 1 : grado_l + 1
                                    for j = 1 : grado_l + 1
                                        for k = 1 : grado_l + 1
                                            w = obj.elemento.ww(k);
                                            Ts(i,j) = Ts(i,j) + w * ss_hg * Lag_eval(i,k) * Lag_eval(j,k);
                                        end
                                    end
                                end
                                Ts = Ts * Jac(e);
                                A((g - 1) * M_nodos + idx, (h - 1) * M_nodos+idx) = ...
                                    A((g - 1) * M_nodos + idx, (h - 1) * M_nodos+idx) - Ts;
                            end
                        end
                    end
                end

                for g = 1 : G
                    for h = 1 : G
                        nf_h = nu_sigma_f(h, mat);
                        if chi(g) ~= 0 && nf_h ~= 0
                            Te = zeros(grado_l + 1);
                            for i = 1 : grado_l + 1
                                for j = 1 : grado_l + 1
                                    for k = 1 : grado_l + 1
                                        w = obj.elemento.ww(k);
                                        Te(i,j) = Te(i,j) + w * chi(g) * nf_h * Lag_eval(i,k) * Lag_eval(j,k);
                                    end
                                end
                            end
                            Te = Te * Jac(e);
                            B((g - 1) * M_nodos + idx, (h - 1) * M_nodos+idx) = ...
                                B((g - 1) * M_nodos + idx, (h - 1) * M_nodos+idx) + Te;
                        end
                    end
                end
            end

            obj.A = A;
            obj.B = B;
        end

        function obj = aplicar_cc(obj)
            G = obj.G;
            M = obj.malla.N * obj.malla.grado_l + 1;
            for g = 1 : G
                base = (g - 1) * M;
                obj.A(base + 1, :) = 0; 
                obj.A(base + 1, base + 1) = 1; 
                obj.A(base + M, :) = 0; 
                obj.A(base + M, base + M) = 1;
                obj.B(base + 1, :) = 0;
                obj.B(base + M, :) = 0;
            end
        end

        function obj = resolver_autovalor(obj, nm)
            [phi, K] = eigs(obj.B, obj.A, nm, 'largestreal');
            obj.keff = diag(K);
            for i = 1:nm
                if i == 1
                    phi(:,i) = abs(phi(:,i));
                end
                phi(:,i) = phi(:,i) / max(phi(:,i));
            end
            obj.phi_inc = phi;
        end
    end
end