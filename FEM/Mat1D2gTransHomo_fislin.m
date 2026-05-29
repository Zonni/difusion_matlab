function params = Mat1D2gTransHomo_fislin()
% Reactor homogeneo 2G sin precursores con perturbacion lineal en fision

    params.nombre = 'Mat1D2gTransHomo_fislin';
    params.G = 2;

    %% Geometria
    params.L = 40;
    params.N = 10;
    params.grado_l = 3;
    params.tam_celdas = params.L / params.N * ones(1, params.N);
    params.materiales = ones(1, params.N);

    %% Propiedades de los materiales
    params.D_0 = [1.32; 0.2772];
    params.sigma_a_0 = [0.0026562; 0.071596];
    params.nu_sigma_f_0 = [0.0074527; 0.08391];

    params.sigma_s_0 = zeros(2,2,1);
    params.sigma_s_0(2,1,1) = 0.023106;

    params.chi = [1.0; 0.0];

    %% Parametros cineticos
    params.beta = [];
    params.lambda = [];
    params.v = [1.27e7; 2.5e5];

    %% Parametros temporales
    params.t_total = 0.5;
    params.pert_D = 0;
    params.pert_a = 0;
    params.pert_f = 1e-7;          
    params.tipo_pert = 'lineal';

    %% Funciones temporales
    params.D_func = cell(2,1);
    params.D_func{1,1} = @(t) params.D_0(1);
    params.D_func{2,1} = @(t) params.D_0(2);

    params.sigma_a_func = cell(2,1);
    params.sigma_a_func{1,1} = @(t) params.sigma_a_0(1);
    params.sigma_a_func{2,1} = @(t) params.sigma_a_0(2);

    params.nu_sigma_f_func = cell(2,1);
    params.nu_sigma_f_func{1,1} = @(t) params.nu_sigma_f_0(1);
    params.nu_sigma_f_func{2,1} = @(t) params.nu_sigma_f_0(2) + params.pert_f * t;

    params.sigma_s_func = cell(2,2,1);
    params.sigma_s_func{2,1,1} = @(t) 0.023106;

    %% Opciones del integrador
    params.opciones_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Stats', 'off');
end