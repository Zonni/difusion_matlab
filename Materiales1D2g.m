% CLASE: Materiales1D2g
% Gestiona la asignación y validación de propiedades físicas para el problema 
% de difusión de neutrones de 2 GRUPOS.
% Las propiedades (D, Sigma) se almacenan como matrices donde:
% - Filas = Grupo de energía (1: Rápido, 2: Térmico).
% - Columnas = Tipo de material.

classdef Materiales1D2g
    properties
        malla           % Referencia a la discretización espacial
        ind_materiales  % Vector (1xN) que asigna un índice de material a cada celda/elemento
        D               % Matriz Coef. Difusión [2 x Tipos_Material]
        sigma_a         % Matriz Sección Eficaz Absorción [2 x Tipos_Material]
        nu_sigma_f      % Matriz Producción de Fisión [2 x Tipos_Material]
        sigma_s12       % Vector Scattering G1->G2 [1 x Tipos_Material]
    end
    
    methods
        function s = Materiales1D2g(materiales, D, sigma_a, nu_sigma_f, sigma_s12, malla)
            
            s.malla = malla;
            N = s.malla.N;
            
            % 1. Validación de la distribución espacial de materiales
            if length(materiales) ~= N
                error('Tamaño inválido: el vector "materiales" debe tener longitud N = %d (número de elementos).', N);
            end
            
            % 2. Obtención de dimensiones de las propiedades
            [nD, mD] = size(D);
            [nA, mA] = size(sigma_a);
            [nF, mF] = size(nu_sigma_f);
            mS = length(sigma_s12);
            
            % 3. Verificación de Grupos de Energía (Deben ser 2 filas)
            if nD ~= 2 || nA ~= 2 || nF ~= 2
                error('Tamaño inválido: las matrices D, sigma_a y nu_sigma_f deben tener dos filas (Grupo 1 y Grupo 2).');
            end
            
            % 4. Verificación de Consistencia de Tipos de Materiales (Columnas iguales)
            if mD ~= mA || mD ~= mF || mD ~= mS
                error('Tamaño inválido: el número de materiales (columnas) difiere entre las propiedades D, sigma_a, nu_sigma_f y sigma_s12.');
            end
            
            % 5. Asignación final de propiedades
            s.ind_materiales = materiales;
            s.D = D;
            s.sigma_a = sigma_a;
            s.nu_sigma_f = nu_sigma_f;
            s.sigma_s12 = sigma_s12;
        end
    end
end