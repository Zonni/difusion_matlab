% CLASE: ElementoFinito
% Define las propiedades del elemento de referencia unidimensional en el intervalo [-1, 1].
% Genera simbólicamente las funciones de forma (Polinomios de Lagrange) y sus derivadas,
% y establece los puntos de integración numérica (Cuadratura de Gauss).

classdef ElementoFinito
    properties
        ref_nodos   % Coordenadas de los nodos en [-1, 1]
        Lag         % Vector de polinomios de Lagrange simbólicos (Funciones de forma Li)
        dLag        % Derivadas simbólicas de los polinomios (dLi/dx)
        xx          % Puntos de integración (Gauss)
        ww          % Pesos de integración (Gauss)
    end
    
    methods
        function s = ElementoFinito(grado_l)
            % 1. Definición de nodos en el elemento de referencia [-1, 1]
            s.ref_nodos = linspace(-1,1,grado_l+1);
            
            % 2. Construcción simbólica de Funciones de Forma (Lagrange)
            syms x
            Lag = sym(zeros(grado_l+1,1));
            
            for i=1:grado_l+1
                Li = sym(1);
                % Productorio clásico de Lagrange: L_i(x) = product((x-x_k)/(x_i-x_k))
                for k=1:grado_l+1
                    if k~=i
                        Li = Li * (x - s.ref_nodos(k)) / (s.ref_nodos(i) - s.ref_nodos(k));
                    end
                end
                Lag(i) = Li;
            end
            s.Lag = Lag;
            
            % 3. Cálculo de derivadas
            s.dLag = diff(Lag, x);
            
            % 4. Obtención de cuadratura de Gauss (Integración numérica exacta)
            % Se solicitan suficientes puntos para integrar el grado del polinomio
            [s.xx, s.ww] = GaussQuad(grado_l+1);
        end
    end
end