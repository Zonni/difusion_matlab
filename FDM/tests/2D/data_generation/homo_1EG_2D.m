%% homo_1EG_2D_generate
% Computes the validated solution for a 2D homogeneous 1-group reactor
% and saves it to "homo_1EG_2D_exact_results.mat" for use in the test suite.
%
% The .mat file stores:
%   k_analytic    — closed-form eigenvalue for a bare square reactor
%   N_array       — mesh sizes used for validation
%   k_eff_num_vec — numerical k-eff at each mesh size (for reference)

clearvars; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'src', '2D'));

%% Problem parameters
L = 350; % cm — full reactor side length
D_val = 0.776; % cm
SigA_val = 0.0244; % 1/cm
NuSigF_val = 0.0260; % 1/cm

% 1-group library arrays (1 material x 1 group)
D_lib = D_val;
sigma_a_lib = SigA_val;
nu_sigf_lib = NuSigF_val;
chi_lib = 1.0;
sigma_s_lib = 0.0;

%% Analytic solution
Bg2 = 2 * (pi / L)^2;
k_analytic = NuSigF_val / (D_val * Bg2 + SigA_val);

fprintf('Analytic k-eff: %.10f\n\n', k_analytic);

%% Numerical solutions at high resolution
N_array = [40, 80, 120, 160];
k_eff_num_vec = zeros(length(N_array), 1);

fprintf('%-12s  %-12s  %-12s  %-12s\n', 'Grid', 'Matrix size', 'k_numeric', 'err (pcm)');
fprintf('%s\n', repmat('-', 1, 56));

for i = 1:length(N_array)
    N = N_array(i);

    % Mesh: N uniform cells per side
    region_lengths_x = (L / N) * ones(N, 1);
    cells_per_region_x = ones(N, 1);

    region_lengths_y = (L / N) * ones(N, 1);
    cells_per_region_y = ones(N, 1);


    mesh  = Mesh_2D_FDM(region_lengths_x, cells_per_region_x, region_lengths_y, cells_per_region_y);

    % Material grid: all cells active, single material
    region_materials = ones(N, N); % (Ny x Nx)

    nuclear_data = NuclearData_2D(region_materials, D_lib, sigma_a_lib, ...
                                  nu_sigf_lib, chi_lib, sigma_s_lib);

    solver = Solver_2D_FDM(mesh, nuclear_data);
    solver = solver.assembleMatrices();
    solver = solver.solveEigenvalues(1);

    k_eff_num_vec(i) = solver.keff(1);
    err_pcm = abs(k_analytic - k_eff_num_vec(i)) * 1e5;

    fprintf('%-12s  %-12d  %-12.10f  %-12.2f\n', ...
        sprintf('%dx%d', N, N), N^2, k_eff_num_vec(i), err_pcm);
end

%% Save results
save_path = fullfile(fileparts(mfilename('fullpath')), '..', 'data', ...
                     'homo_1EG_2D_exact_results.mat');
save(save_path, 'N_array', 'k_eff_num_vec', 'k_analytic');
fprintf('\nResults saved to: %s\n', save_path);