classdef Malla1D
    properties
        L
        N
        grado_l
        x_nodos
        nodos
        tamano_celdas
        x_celdas
    end
    
    methods
        function s = Malla1D(L, N, grado_l, tamano_celdas)
            s.L = L;
            s.N = N;
            s.grado_l = grado_l;

            s.nodos = zeros(N , grado_l+1);
            for e = 1:N
                s.nodos(e,:) = (e-1)*grado_l + (1:grado_l+1);
            end
            
            if sum(tamano_celdas) ~= L || length(tamano_celdas) ~= N
                error('Tamaño Invalido: el vector tamaño debe tener longitud L = %d y sus elementos deben sumar N = %d.', L, N);
            else
                s.tamano_celdas = tamano_celdas;
            end
            
            x_nodos = linspace(0, tamano_celdas(1), grado_l + 1);  
            x_e = tamano_celdas(1);       
            for e = 2:N
                x_local = linspace(0, tamano_celdas(e), grado_l + 1);        
                x_elem = x_e + x_local;   
                
                x_nodos = [x_nodos, x_elem(2:end)];
                x_e = x_e + tamano_celdas(e);
            end
            s.x_nodos = x_nodos;

            s.x_celdas=[0, cumsum(tamano_celdas)];
        end
    end
end
