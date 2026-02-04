% CLASE: Materiales1D1g
% Gestiona las propiedades físicas para el problema de 1 GRUPO de energía.
% Asigna a cada celda de la malla un material específico y almacena sus
% secciones eficaces (D, Sigma_a, Nu_Sigma_f).

classdef Materiales1D1g
    properties
        malla           % Referencia a la discretización espacial
        ind_materiales  % Vector (1xN): Indica el Ind del material presente en cada celda
        D               % Vector Coef. Difusión [1 x Tipos_Material]
        sigma_a         % Vector Sección Absorción [1 x Tipos_Material]
        nu_sigma_f      % Vector Producción Fisión [1 x Tipos_Material]
    end
    
    methods
        function s = Materiales1D1g(materiales, D, sigma_a, nu_sigma_f, malla)
            s.malla = malla;
            N = s.malla.N;
            
            % Validación: El vector 'materiales' debe definir el tipo de material para cada una de las N celdas
            if  length(materiales) ~= N 
                error('Tamaño Inválido: el vector materiales debe tener longitud N = %d.', N);
            else
                % Almacenamiento de propiedades
                s.ind_materiales = materiales;
                s.D = D;
                s.sigma_a = sigma_a;
                s.nu_sigma_f = nu_sigma_f;
            end
        end
        
    end
end