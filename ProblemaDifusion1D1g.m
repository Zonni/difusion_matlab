% CLASE: ProblemaDifusion1D1g
% Esta clase resuelve el problema de Difusión de Neutrones en 1D con 1 grupo de energía.
% Permite la resolución numérica mediante dos métodos: Elementos Finitos (FEM) o 
% Diferencias Finitas (FDM). Su objetivo principal es ensamblar las matrices del sistema,
% aplicar condiciones de contorno y resolver el problema de autovalores para encontrar
% la constante de multiplicación efectiva (keff) y el flujo de neutrones (phi).

classdef ProblemaDifusion1D1g
    properties
        malla           % Objeto con información de la discretización espacial
        materiales      % Propiedades físicas (D, sigma_a, nu_sigma_f)
        elemento        % Definición del elemento finito (funciones de forma, cuadratura)
        A               % Matriz de rigidez/pérdidas (Fugas + Absorción)
        B               % Matriz de producción (Fisión)
        phi_inc         % Vectores propios (Flujos) calculados
        keff            % Autovalores (k efectiva) calculados
    end
    
    methods
        % CONSTRUCTOR: Inicializa el objeto con la malla, materiales y tipo de elemento.
        function s = ProblemaDifusion1D1g(malla, materiales, elemento)
            s.malla = malla;
            s.materiales = materiales;
            s.elemento = elemento;
        end
        
        % MÉTODO FEM: Ensambla las matrices A y B usando el Método de Elementos Finitos.
        function s = ensamblar_matrices_fem(s)
            % 1. Extracción de parámetros de la malla y materiales
            N = s.malla.N;
            grado_l = s.malla.grado_l;
            tamano_celdas = s.malla.tamano_celdas;
            D = s.materiales.D;
            sigma_a = s.materiales.sigma_a;
            nu_sigma_f = s.materiales.nu_sigma_f;
            ind_materiales = s.materiales.ind_materiales;
            
            % 2. Evaluar polinomios de Lagrange y derivadas en puntos de integración
            Lag_eval = double(subs(s.elemento.Lag, sym('x'), s.elemento.xx'));
            dLag_eval = double(subs(s.elemento.dLag, sym('x'), s.elemento.xx'));
            Jac = tamano_celdas/2; % Jacobiano de la transformación
            
            % 3. Matrices globales dispersas
            A = sparse(N*grado_l+1, N*grado_l+1);
            B = sparse(N*grado_l+1, N*grado_l+1);
            
            % 4. Bucle sobre elementos para integración numérica y ensamblaje
            for e = 1:N
                idx = s.malla.nodos(e, :); % Índices globales de los nodos del elemento
                Ke = zeros(grado_l+1);    
                Me = zeros(grado_l+1);     
                Te = zeros(grado_l+1);     
                
                % Integración numérica (Cuadratura de Gauss)
                for i=1:grado_l+1
                    for j=1:grado_l+1
                        for k=1:grado_l+1
                            
                            Ke(i,j) = Ke(i,j) + s.elemento.ww(k)*D(ind_materiales(e))*dLag_eval(i,k)*dLag_eval(j,k);
                           
                            Me(i,j) = Me(i,j) + s.elemento.ww(k)*sigma_a(ind_materiales(e))*Lag_eval(i,k)*Lag_eval(j,k);
                        
                            Te(i,j) = Te(i,j) + s.elemento.ww(k)*nu_sigma_f(ind_materiales(e))*Lag_eval(i,k)*Lag_eval(j,k);
                        end
                    end
                end
                
                % Aplicar Jacobiano a las integrales
                Ke = Ke / Jac(e);
                Me = Me * Jac(e);
                Te = Te * Jac(e);
                
                % Ensamblaje en matrices globales
                A(idx,idx) = A(idx,idx) + Ke + Me;
                B(idx,idx) = B(idx,idx) + Te;
            end
            
            % 5. Guardar matrices en el objeto
            s.A = A;
            s.B = B;
        end
        
        % MÉTODO FDM: Ensambla las matrices A y B usando Diferencias Finitas centradas.
        function s = ensamblar_matrices_fdm(s)
            % 1. Extracción de parámetros
            N = s.malla.N;
            L = s.malla.L;
            D = s.materiales.D;
            sigma_a = s.materiales.sigma_a;
            nu_sigma_f = s.materiales.nu_sigma_f;
            ind_materiales = s.materiales.ind_materiales;
            
            % 2. Definición espacial auxiliar
            x_vec = linspace(0, L, N+2);
            phi_analytic = sin(pi * x_vec/L);
            
            % 3. Reserva de memoria para matrices dispersas
            A = spalloc(N, N, N*3); 
            B = spalloc(N, N, N);  
            
            delta_x = L / (N+1); % Longitud de un elemento de la malla
            
            % 4. Definición de ecuaciones para el primer nodo 
            A(1, 1) = 2*D(ind_materiales(1))/delta_x^2 + sigma_a(ind_materiales(1));
            A(1, 2) = - D(ind_materiales(1))/delta_x^2;
            B(1, 1) = nu_sigma_f(ind_materiales(1));
            
            % 5. Bucle para nodos internos
            for i = 2:N-1
                A(i, i-1) = - D(ind_materiales(i))/delta_x^2;
                A(i, i)   = 2*D(ind_materiales(i))/delta_x^2 + sigma_a(ind_materiales(i));
                A(i, i+1) = - D(ind_materiales(i))/delta_x^2;
                
                B(i, i) = nu_sigma_f(ind_materiales(i));
            end
            
            % 6. Definición de ecuaciones para el último nodo 
            A(N, N-1) = - D(ind_materiales(N))/delta_x^2 ;
            A(N, N)   = 2*D(ind_materiales(N))/delta_x^2 + sigma_a(ind_materiales(N));
            B(N, N)   = nu_sigma_f(ind_materiales(N));
            
            s.A = A;
            s.B = B;
        end

        % MÉTODO CC: Aplica condiciones de contorno de Dirichlet homogéneas (Phi=0 en extremos).
        function s = aplicar_cc(s)
            % Modifica primera fila para forzar Phi_1 = 0
            s.A(1,:) = 0; s.A(1,1)=1; s.B(1,:)=0;
            % Modifica última fila para forzar Phi_end = 0
            s.A(end,:) = 0; s.A(end,end)=1; s.B(end,:)=0;
        end
        
        % MÉTODO SOLVER: Resuelve el problema de autovalores generalizado.
        function s = resolver_autovalor(s,nm)
            % 1. Cálculo de autovalores y autovectores (nm = número de modos)
            [phi,K] = eigs(s.B, s.A, nm);
            s.keff = diag(K);
            
            % 2. Normalización de los flujos
            for i = 1:nm
                % Asegurar signo positivo global
                if phi(2,i)<0
                    phi(:,i)=-phi(:,i);
                end
                % Normalizar al valor máximo unitario
                phi(:,i)=phi(:,i)/max(phi(:,i));
                
            end
            s.phi_inc = phi;
        end

        % MÉTODO GRAFICAR: Visualiza los modos del flujo calculados.
        function s = graficar(s,modos)
            x  = s.malla.x_nodos;
            figure; hold on; grid on;
            estilos_modo = {'-','--',':','-.'}; 
        
            % Bucle para graficar cada modo solicitado
            for idx = 1:length(modos)
                modo = modos(idx);
                estilo = estilos_modo{mod(idx-1,length(estilos_modo))+1};
        
                plot(x, s.phi_inc(1:end,modo), 'b','LineStyle', estilo, 'LineWidth', 1.5);
            end
        
            xlabel('x (cm)'); ylabel('\phi(x)');
            title('Modos FEM');
            
            % Generación dinámica de la leyenda
            leyenda = {};
            for idx = 1:length(modos)
                leyenda{end+1} = sprintf('Modo %d', modos(idx));
            end
            legend(leyenda,'Location','best');
        end
    end
end


