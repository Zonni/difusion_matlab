function params = Mat1D7gTransCte()
% Reactor homogeneo 7G sin precursores sin perturbacion

    params.nombre = 'Mat1D7gTransCte';
    params.G = 7;

    %% Geometria
    params.L = 40;                       
    params.N = 10;
    params.grado_l = 3;
    params.tam_celdas = params.L / params.N * ones(1, params.N);
    params.materiales = ones(1, params.N);   

    %% Propiedades de los materiales (G x 1)
    % Coeficientes de difusion (cm)
    params.D_0 = [1.500e+0;   % G1
                  1.200e+0;   % G2
                  9.000e-1;   % G3
                  6.000e-1;   % G4
                  4.000e-1;   % G5
                  3.000e-1;   % G6
                  2.000e-1];  % G7

    % Secciones eficaces de absorcion (cm^{-1})
    params.sigma_a_0 = [2.500e-3;   % G1
                        3.000e-3;   % G2
                        4.000e-3;   % G3
                        6.000e-3;   % G4
                        1.000e-2;   % G5
                        2.000e-2;   % G6
                        4.000e-2];  % G7

    % nu * sigma_f (cm^{-1})
    params.nu_sigma_f_0 = [2.000e-3;   % G1
                           1.500e-3;   % G2
                           1.000e-3;   % G3
                           5.000e-4;   % G4
                           2.000e-4;   % G5
                           1.000e-4;   % G6
                           0.0];       % G7

    % Scattering de grupo h a grupo g (g,h,mat)
    params.sigma_s_0 = zeros(7, 7, 1); 
    params.sigma_s_0(2, 1, 1) = 0.025; % G1 -> G2
    params.sigma_s_0(3, 2, 1) = 0.020; % G2 -> G3
    params.sigma_s_0(4, 3, 1) = 0.015; % G3 -> G4
    params.sigma_s_0(5, 4, 1) = 0.010; % G4 -> G5
    params.sigma_s_0(6, 5, 1) = 0.008; % G5 -> G6
    params.sigma_s_0(7, 6, 1) = 0.005;  % G6 -> G7

    % Chi
    params.chi = [0.6; 0.3; 0.1; 0.0; 0.0; 0.0; 0.0];

    %% Precursores
    params.beta = [];                
    params.lambda = [];
    params.v = [3.0e7; 2.0e7; 1.0e7; 5.0e6; 2.0e6; 5.0e5; 2.0e5];

    %% Parametros temporales
    params.t_total = 0.5;
    params.pert_D   = 0;
    params.pert_a   = 0;
    params.pert_f   = 0;
    params.tipo_pert = 'ninguna';

    %% Funciones dependientes del tiempo (cell G x 1)
    % Difusion
    params.D_func = cell(7, 1);
    for g = 1 : 7
        params.D_func{g,1} = @(t) params.D_0(g);
    end

    % Absorcion
    params.sigma_a_func = cell(7, 1);
    for g = 1 : 7
        params.sigma_a_func{g,1} = @(t) params.sigma_a_0(g);
    end

    % Fision
    params.nu_sigma_f_func = cell(7, 1);
    for g = 1 : 7
        params.nu_sigma_f_func{g,1} = @(t) params.nu_sigma_f_0(g);
    end

    % Scattering (cell 7 x 7 x 1)
    params.sigma_s_func = cell(7, 7, 1);
    for g = 1 : 7
        for h = 1 : 7
            params.sigma_s_func{g,h,1} = [];
        end
    end
    params.sigma_s_func{2,1,1} = @(t) 0.025;
    params.sigma_s_func{3,2,1} = @(t) 0.020;
    params.sigma_s_func{4,3,1} = @(t) 0.015;
    params.sigma_s_func{5,4,1} = @(t) 0.010;
    params.sigma_s_func{6,5,1} = @(t) 0.008;
    params.sigma_s_func{7,6,1} = @(t) 0.005;

    %% Opciones del integrador ODE
    params.opciones_ode = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Stats', 'off');
end