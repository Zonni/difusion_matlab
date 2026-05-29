clearvars;
close all;

%% MATERIALS
num_materials = 2;
num_groups    = 2;

%            Group:  [  fast,  thermal ]
D_lib         = [1.440,  0.390;   % Material 1: Fuel
                  1.720,  0.250];  % Material 2: Reflector

sigma_a_lib   = [0.0095, 0.0800;  % Fuel
                  0.0002, 0.0100]; % Reflector

nu_sigma_f_lib = [0.000,  0.1325;  % Fuel: fission only in thermal group
                   0.000,  0.000];  % Reflector: no fission

chi_lib        = [1, 0;            % Fuel: all neutrons born fast
                   1, 0];           % Reflector: (irrelevant, no fission)

% sigma_s_lib(material, g_from, g_to)
% Only fast->thermal downscattering; no upscattering
sigma_s_lib = zeros(num_materials, num_groups, num_groups);
sigma_s_lib(1, 1, 2) = 0.0180;   % Fuel:      fast -> thermal
sigma_s_lib(2, 1, 2) = 0.0430;   % Reflector: fast -> thermal

% Region-to-material map: Reflector | Fuel | Reflector
region_materials = [2; 1; 2];

%% GEOMETRY
region_lengths   = [30;  80;  30];  % cm
cells_per_region = [20;  50;  20];

%% SOLVE
mesh      = Mesh_1D_FDM(region_lengths, cells_per_region);
materials = Materials(region_materials, D_lib, sigma_a_lib, ...
                      nu_sigma_f_lib, chi_lib, sigma_s_lib);

solver = Solver_1D_FDM(mesh, materials);
solver = solver.assembleMatrices();
solver = solver.solveEigenvalues(3);
%solver.displayProblem();
solver.plotPhi();