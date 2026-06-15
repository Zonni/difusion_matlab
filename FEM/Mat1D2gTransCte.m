function params = Mat1D2gTransCte()
% Reactor homogeneo 2G sin precursores sin perturbacion

    params.nombre = 'Mat1D2gTransCte';
    params.G = 2;                     

    %% Geometria
    params.L = 40;
    params.N = 10;
    params.grado_l = 3;
    params.tam_celdas = params.L / params.N * ones(1, params.N);
    params.materiales = ones(1, params.N);

    %% Propiedades de los materiales (G x n_mat)

    % Coeficientes de difusion
    params.D_0 = [1.32;       % grupo 1 (rapido)
                  0.2772];    % grupo 2 (termico)

    % Secciones eficaces de absorcion
    params.sigma_a_0 = [0.0026562;   % grupo 1
                        0.071596];   % grupo 2

    % Secciones eficaces de produccion (nu * sigma_f)
    params.nu_sigma_f_0 = [0.0074527;   % grupo 1
                           0.08391];    % grupo 2

    % Matriz de scattering  sigma_s(g, h, material) de grupo h a grupo g
    % Dimension: G x G x n_mat
    params.sigma_s_0 = zeros(2, 2, 1);
    params.sigma_s_0(2, 1, 1) = 0.023106;   

    % Chi
    params.chi = [1.0; 0.0];        % G x 1

    %% Parametros cineticos
    params.beta   = [];              
    params.lambda = [];
    params.v = [1.27e7; 2.5e5];      % Velocidades por grupo (cm/s)

    %% Parametros temporales
    params.t_total = 0.5;
    params.pert_D   = 0;
    params.pert_a   = 0;
    params.pert_f   = 0;
    params.tipo_pert = 'ninguna';

    %% Funciones temporales (cell arrays G x n_mat)

    % Difusion
    params.D_func = cell(2, 1);
    params.D_func{1,1} = @(t) params.D_0(1);
    params.D_func{2,1} = @(t) params.D_0(2);

    % Absorcion
    params.sigma_a_func = cell(2, 1);
    params.sigma_a_func{1,1} = @(t) params.sigma_a_0(1);
    params.sigma_a_func{2,1} = @(t) params.sigma_a_0(2);

    % Produccion (nu*sigma_f)
    params.nu_sigma_f_func = cell(2, 1);
    params.nu_sigma_f_func{1,1} = @(t) params.nu_sigma_f_0(1);
    params.nu_sigma_f_func{2,1} = @(t) params.nu_sigma_f_0(2);

    % Scattering (cell G x G x n_mat)
    params.sigma_s_func = cell(2, 2, 1);
    for g = 1:2
        for h = 1:2
            params.sigma_s_func{g, h, 1} = [];
        end
    end
    params.sigma_s_func{2, 1, 1} = @(t) 0.023106;

    %% Opciones del integrador
    params.opciones_ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Stats', 'off');
end