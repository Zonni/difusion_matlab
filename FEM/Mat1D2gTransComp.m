function params = Mat1D2gTransComp()
% Reactor heterogeneo 2G con precursores y perturbacion lineal en absorcion y fision

    params.nombre = 'Mat1D2gTransComp';
    params.G = 2;

    %% Geometria
    params.L = 350;                         
    params.N = 14;                           
    params.grado_l = 3;
    params.tam_celdas = [5, 25 * ones(1,12), 45];
    params.materiales = [1, 2 * ones(1,12), 1]; 

    %% Propiedades de los materiales (G x n_mat)
    % Material 1 (reflector)
    D1_r = 1.5;        D2_r = 0.4;
    sa1_r = 0.01;      sa2_r = 0.085;    s12_r = 0.02;
    nf1_r = 0.0;       nf2_r = 0.0;

    % Material 2 (combustible)
    D1_f = 1.32;       D2_f = 0.2772;
    sa1_f = 0.0026562; sa2_f = 0.071596; s12_f = 0.023106;
    nf1_f = 0.0074527; nf2_f = 0.08391;

    params.D_0 = [D1_r, D1_f; D2_r, D2_f];
    params.sigma_a_0 = [sa1_r, sa1_f; sa2_r, sa2_f];
    params.nu_sigma_f_0 = [nf1_r, nf1_f; nf2_r, nf2_f];

    % Scattering para cada material
    params.sigma_s_0 = zeros(2, 2, 2);
    params.sigma_s_0(2, 1, 1) = s12_r;   % material 1
    params.sigma_s_0(2, 1, 2) = s12_f;   % material 2

    params.chi = [1.0; 0.0];

    %% Parametros cineticos (con precursores)
    params.beta   = [0.000247, 0.0013845, 0.001222, 0.0026455, 0.000832, 0.000169];
    params.lambda = [0.0127,   0.0317,    0.115,    0.311,     1.4,     3.87];
    params.v = [1.27e7; 2.5e5];

    %% Parametros temporales
    params.t_total = 0.5;          
    params.pert_D = 0;
    params.pert_a = 1e-3;          
    params.pert_f = 1e-3;          
    params.tipo_pert = 'lineal';

    %% Funciones temporales (cell G x n_mat)
    
    params.D_func = cell(2,2);
    params.D_func{1,1} = @(t) D1_r;
    params.D_func{2,1} = @(t) D2_r;
    params.D_func{1,2} = @(t) D1_f;
    params.D_func{2,2} = @(t) D2_f;

    params.sigma_a_func = cell(2,2);
    params.sigma_a_func{1,1} = @(t) sa1_r;
    params.sigma_a_func{2,1} = @(t) sa2_r;
    params.sigma_a_func{1,2} = @(t) sa1_f;
    params.sigma_a_func{2,2} = @(t) sa2_f + params.pert_a;

    params.nu_sigma_f_func = cell(2,2);
    params.nu_sigma_f_func{1,1} = @(t) nf1_r;
    params.nu_sigma_f_func{2,1} = @(t) nf2_r;
    params.nu_sigma_f_func{1,2} = @(t) nf1_f;
    params.nu_sigma_f_func{2,2} = @(t) nf2_f + params.pert_f; 

    % Scattering
    params.sigma_s_func = cell(2,2,2);
    % Material 1
    params.sigma_s_func{2,1,1} = @(t) s12_r;   
    % Material 2
    params.sigma_s_func{2,1,2} = @(t) s12_f;   

    %% Opciones del integrador ODE
    params.opciones_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Stats', 'off');
end