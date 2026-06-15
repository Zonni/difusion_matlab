function params = Mat1D1gTransHomo_fisesc()
% Reactor homogeneo 1G sin precursores, perturbacion escalon en fision

    params.nombre = 'Mat1D1gTransHomo_fisesc';
    params.G = 1;

    %% Geometria
    params.L = 350;
    params.N = 14;
    params.grado_l = 3;
    params.tam_celdas = 25 * ones(1, 14);
    params.materiales = ones(1, 14);

    %% Propiedades
    params.D_0 = 0.776;
    params.sigma_a_0 = 0.0244;
    params.nu_sigma_f_0 = 0.0260;
    params.sigma_s_0 = zeros(1, 1, 1);
    params.chi = 1.0;

    %% Cinetica
    params.beta = [];
    params.lambda = [];
    params.v = 1.25e7;

    %% Transitorio
    params.t_total = 0.5;
    params.pert_D = 0;
    params.pert_a = 0;
    params.pert_f = 1e-7;                  
    params.tipo_pert = 'escalon';

    %% Funciones
    params.D_func = cell(1, 1);
    params.D_func{1,1} = @(t) params.D_0;

    params.sigma_a_func = cell(1, 1);
    params.sigma_a_func{1,1} = @(t) params.sigma_a_0;

    params.nu_sigma_f_func = cell(1, 1);
    params.nu_sigma_f_func{1,1} = @(t) params.nu_sigma_f_0 + params.pert_f;
    
    % Scattering
    params.sigma_s_func = cell(1, 1, 1);

    %% Opciones
    params.opciones_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Stats', 'off');
end