function check_banco(banco_datos, tol, bool_graf, metodo)
% Comprueba la consistencia del metodo numerico comparando con solucion analitica

    for i = 1 : length(banco_datos)
        banco = banco_datos{i};
        
        % Cargar parametros del banco actual
        params = cargar_params(banco);
        
        % Ejecutar simulacion segun el metodo elegido
        switch metodo
            case 'fem_gTrans'
                [transitorio, t, ~] = fem_gTrans(params, bool_graf);
            case 'fem_gTransPC'
                [transitorio, t, ~] = fem_gTransPC(params, bool_graf);
            case 'fem_2gTransPC'
                if params.G ~= 2
                    error('El metodo fem_2gTransPC solo admite 2 grupos de energia');
                end
                [transitorio, t, ~] = fem_2gTransPC(params, bool_graf);
            otherwise
                error('Metodo no reconocido: %s', metodo);
        end
        Pot_num = transitorio.Pot_num;
        
        % Determinar si el caso admite solucion analitica
        es_homogeneo = (isscalar(unique(params.materiales)));
        sin_precursores = isempty(params.beta);
        
        % Casos analizables segun el numero de grupos
        if params.G == 1
            % En 1G todos los tipos de perturbacion constantes, escalon y lineales tienen solucion analitica
            caso_analizable = es_homogeneo && sin_precursores && ...
                              ismember(params.tipo_pert, {'ninguna', 'escalon', 'lineal'});
        elseif params.G == 2
            % En 2G todos los tipos de perturbacion constantes y escalon tienen solucion analitica
            caso_analizable = es_homogeneo && sin_precursores && ...
                              ismember(params.tipo_pert, {'ninguna', 'escalon'});
        else
            caso_analizable = false;   % No hay soluciones para G > 2
        end

        if ~caso_analizable
            fprintf('-> Sin solucion analitica\n');
        else
            % Calcular solucion analitica segun el numero de grupos y tipo de perturbacion
            if params.G == 1
                Pot_an = pot_an_1g(params, t);
                if strcmp(params.tipo_pert, 'ninguna')
                    error = abs(Pot_num(end) - 1);
                else
                    error = abs( Pot_num(end) - Pot_an(end) );
                end
            elseif params.G == 2
                if strcmp(params.tipo_pert, 'ninguna')
                    Pot_an = ones(size(t));
                    error = abs(Pot_num(end) - 1);
                elseif strcmp(params.tipo_pert, 'escalon')
                    Pot_an = pot_an_2g(params, t);
                    error = abs( Pot_num(end) - Pot_an(end) );
                else
                    error('Tipo de perturbacion desconocido para 2G: %s', params.tipo_pert);
                end
            else
                error('Grupos no soportados en solucion analitica: %d', params.G);
            end
            
            % Mostrar resultados y verificar tolerancia
            assert(error < tol, ...
                    '-> Banco %s: error = %.3e > %.3e', banco, error, tol);
            fprintf('-> Banco %s: OK (error = %.3e)\n', banco, error);
            
            % Graficas comparativas (opcional: saveas)
            if bool_graf
                figure('Name', sprintf('Comparacion %s', banco));
                plot(t, Pot_num, 'b-', 'LineWidth', 2);
                hold on;
                plot(t, Pot_an, 'ro', 'LineWidth', 2);
                xlabel('Tiempo (s)');
                ylabel('Potencia normalizada');
                legend('FEM', 'Analitica', 'Location', 'best');
                grid on;
                hold off;

                saveas(gcf, ['Comparacion_', banco, '.eps'], 'epsc');
            end
        end
    end 
end

% -------------------------------------------------------------------------
function Pot_an = pot_an_1g(params, t)
% Solucion analitica para un grupo de energia
% Soporta perturbaciones constantes, escalon y lineales en absorcion y fision

    % Buckling geometrico
    B2 = (pi / params.L)^2;
    % Secciones eficaces iniciales
    D0   = params.D_0(1);
    Sa0  = params.sigma_a_0(1);
    nf0  = params.nu_sigma_f_0(1);
    v    = params.v(1);

    % Constante de multiplicacion efectiva inicial (analitica)
    keff0_an = nf0 / (D0 * B2 + Sa0);

    switch params.tipo_pert
        case 'ninguna'
            Pot_an = ones(size(t));

        case 'escalon'
            if params.pert_a ~= 0
                % Perturbacion escalon en absorcion
                Pot_an = exp( - v * params.pert_a * t );
            elseif params.pert_f ~= 0
                % Perturbacion escalon en fision
                Pot_an = exp( (v * params.pert_f / keff0_an) * t );
            else
                Pot_an = ones(size(t));
            end

        case 'lineal'
            if params.pert_a ~= 0
                % Perturbacion lineal en absorcion
                Pot_an = exp( - (v * params.pert_a / 2) * t.^2 );
            elseif params.pert_f ~= 0
                % Perturbacion lineal en fision
                Pot_an = exp( (v * params.pert_f / (2 * keff0_an)) * t.^2 );
            else
                Pot_an = ones(size(t));
            end

        otherwise
            error('Tipo de perturbacion no soportado en 1G: %s', params.tipo_pert);
    end
end

% -------------------------------------------------------------------------
function Pot_an = pot_an_2g(params, t)
% Solucion analitica para dos grupos de energia
% Soporta perturbaciones constantes y escalon en absorcion y fision

    % Constantes
    B = pi / params.L;
    v1 = params.v(1);   v2 = params.v(2);
    D1 = params.D_0(1); D2 = params.D_0(2);
    Sa1 = params.sigma_a_0(1);
    Sa2_0 = params.sigma_a_0(2);
    S12 = params.sigma_s_0(2, 1, 1);
    nusF1 = params.nu_sigma_f_0(1);
    nusF2_0 = params.nu_sigma_f_0(2);

    % Determinar tipo de perturbacion
    if params.pert_f ~= 0
        delta_nusF2 = params.pert_f;
        nusF2 = nusF2_0 + delta_nusF2;
        Sa2 = Sa2_0;
        delta_Sa2 = 0;
        tipo_pert = 'fision';
    elseif params.pert_a ~= 0
        delta_Sa2 = params.pert_a;
        nusF2 = nusF2_0;
        Sa2 = Sa2_0 + delta_Sa2;
        delta_nusF2 = 0;
        tipo_pert = 'absorcion';
    else
        nusF2 = nusF2_0;
        Sa2 = Sa2_0;
        delta_nusF2 = 0;
        delta_Sa2 = 0;
        tipo_pert = 'ninguna';
    end
    
    % Coeficientes
    C = (D2 * B^2 + Sa2_0) / S12;
    switch tipo_pert
        case 'fision'
            M = v1 * (D1 * B^2 + Sa1 + S12 - nusF1) + v2 * (D2 * B^2 + Sa2_0);
            N = - v1 * v2 * S12 * delta_nusF2;
        case 'absorcion'
            M = v2 * (D2 * B^2 + Sa2_0 + delta_Sa2) + v1 * (nusF2_0 / C);
            N = v1 * v2 * (nusF2_0 / C) * delta_Sa2;
        otherwise
            M = 0; 
            N = 0;
    end

    disc = M^2 - 4*N;
    beta1 = (-M + sqrt(disc)) / 2;
    beta2 = (-M - sqrt(disc)) / 2;

    T1 = zeros(size(t));
    T2 = zeros(size(t));

    if strcmp(tipo_pert, 'fision')
        denom = beta2 - beta1;
        D0 = D2 * B^2 + Sa2_0;
        for i = 1 : length(t)
            e1 = exp(beta1 * t(i));
            e2 = exp(beta2 * t(i));
            T2(i) = (beta2*e1 - beta1*e2) / denom;
            T1(i) = ( (beta1 * beta2 / (v2 * denom))*(e1 - e2) + (D0 / denom)*(beta2 * e1 - beta1 * e2) ) / D0;
        end
    elseif strcmp(tipo_pert, 'absorcion')
        denom = beta2 - beta1;
        v2_dSa = v2 * delta_Sa2;
        D0 = D2*B^2 + Sa2_0;
        Dpert = D0 + delta_Sa2;
        for i = 1 : length(t)
            e1 = exp(beta1 * t(i));
            e2 = exp(beta2 * t(i));
            T2(i) = ((beta2 + v2_dSa)*e1 - (beta1 + v2_dSa)*e2) / denom;
           
            num1 = (beta2 + v2_dSa) * (beta1 + v2 * Dpert) * e1;
            num2 = (beta1 + v2_dSa) * (beta2 + v2 * Dpert) * e2;
            T1(i) = (num1 - num2) / (v2 * D0 * denom);
        end
    else
        T1(:) = 1;
        T2(:) = 1;
    end

    num   = nusF1 * T1 + nusF2 * T2;
    dem = nusF1 * 1  + nusF2 * 1;
    Pot_an = num / dem;
end