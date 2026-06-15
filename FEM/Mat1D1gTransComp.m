function params = Mat1D1gTransComp()
% Reactor heterogeneo 1G con precursores y perturbaciones lineales

    params.nombre = 'Mat1D1gTransComp';
    params.G = 1;

    %% Geometria
    params.L = 350;
    params.N = 14;
    params.grado_l = 3;
    params.tam_celdas = [5, 25 * ones(1, 12), 45];
    params.materiales = [1, 2 * ones(1, 12), 1];   % 2 materiales

    %% Propiedades de los materiales (G x n_mat, n_mat=2)
    % Material 1 (reflector)
    D_r = 1.446;   sa_r = 0.0077;   nf_r = 0.0;
    % Material 2 (combustible)
    D_f = 0.776;   sa_f = 0.0244;   nf_f = 0.0260;

    params.D_0 = [D_r, D_f];
    params.sigma_a_0 = [sa_r, sa_f];
    params.nu_sigma_f_0 = [nf_r, nf_f];

    % Scattering
    params.sigma_s_0 = zeros(1, 1, 2);
    % Chi
    params.chi = 1.0;

    %% Parametros cineticos (con precursores)
    params.beta   = [0.000247, 0.0013845, 0.001222, 0.0026455, 0.000832, 0.000169];
    params.lambda = [0.0127,   0.0317,    0.115,    0.311,     1.4,     3.87];
    params.v = 1.25e7;

    %% Parametros temporales
    params.t_total = 0.5;
    params.pert_D = 0;
    params.pert_a = 1e-3;
    params.pert_f = 1e-3;
    params.tipo_pert = 'lineal';          

    %% Funciones temporales (cell G x n_mat)
    % Material 1 (constante)
    params.D_func = cell(1, 2);
    params.D_func{1,1} = @(t) D_r;
    params.D_func{1,2} = @(t) D_f + params.pert_D * t;

    params.sigma_a_func = cell(1, 2);
    params.sigma_a_func{1,1} = @(t) sa_r;
    params.sigma_a_func{1,2} = @(t) sa_f + params.pert_a * t;

    params.nu_sigma_f_func = cell(1, 2);
    params.nu_sigma_f_func{1,1} = @(t) nf_r;
    params.nu_sigma_f_func{1,2} = @(t) nf_f + params.pert_f * t;

    % Scattering
    params.sigma_s_func = cell(1, 1, 2);

    %% Opciones
    params.opciones_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Stats', 'off');
end