%% homo_2EG_2D_generate
% Computes the validated solution for a 2D homogeneous 2-group reactor
% and saves it to "homo_2EG_2D_exact_results.mat" for use in the test suite.
%
% Analytical solution is derived from the 3D lambda-modes formulation
% (Appendix A) by dropping the z-dimension (2D square reactor):
%
%   B^2_{m,n} = (m*pi/Lx)^2 + (n*pi/Ly)^2
%
%   psi_2(x,y) = k * sin(m*pi*x/Lx) * sin(n*pi*y/Ly)
%   psi_1(x,y) = (D2*B^2 + sigma_a2) / sigma_s12 * psi_2(x,y)
%
%   lambda_{m,n} = [nu_sigf1*(D2*B^2 + sigma_a2) + nu_sigf2*sigma_s12]
%               / [(D2*B^2 + sigma_a2)*(sigma_a1 + sigma_s12 + D1*B^2)]
%
%   Normalisation (2D, integral of |sin*sin| over [0,Lx]x[0,Ly] = 4*Lx*Ly/pi^2
%   for m=n=1):
%
%   k = (pi^2 / 4) * sigma_s12
%       / [sigma_f1*(D2*B^2 + sigma_a2) + sigma_s12*sigma_f2]
%
% The .mat file stores:
%   k_analytic       — analytic fundamental eigenvalue (m=n=1)
%   phi1_func        — function handle: fast flux  psi_1(x,y)
%   phi2_func        — function handle: thermal flux psi_2(x,y)
%   N_array          — mesh sizes used for the convergence study
%   k_eff_num_vec    — numerical k-eff at each mesh size

clearvars; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'src', '2D'));

%% ========================================================================
%  1. PROBLEM PARAMETERS — 2-group homogeneous bare square reactor
% =========================================================================

L = 200;   % cm — full reactor side length (Lx = Ly = L)

% --- Nuclear data (2 materials x 2 groups) --------------------------------
%                  group 1 (fast)   group 2 (thermal)
D_lib        = [   1.500,            0.400  ];   % (1 x 2)  diffusion coeff (cm)
sigma_a_lib  = [   0.010,            0.080  ];   % (1 x 2)  absorption xs (1/cm)
nu_sigf_lib  = [   0.000,            0.185  ];   % (1 x 2)  nu*fission xs (1/cm)
sigma_f_lib  = [   0.000,            0.065  ];   % (1 x 2)  fission xs (1/cm) — for normalisation only
chi_lib      = [   1.000,            0.000  ];   % (1 x 2)  fission spectrum

% Scattering matrix (1 material x 2 groups x 2 groups)
% sigma_s_lib(mat, g_from, g_to)
% Only downscatter 1->2 is non-zero
sigma_s_lib        = zeros(1, 2, 2);
sigma_s_lib(1,1,2) = 0.025;   % fast -> thermal (1/cm)

% Shorthand scalars for analytical formulae
D1       = D_lib(1);
D2       = D_lib(2);
siga1    = sigma_a_lib(1);
siga2    = sigma_a_lib(2);
nusigf1  = nu_sigf_lib(1);
nusigf2  = nu_sigf_lib(2);
sigf1    = sigma_f_lib(1);
sigf2    = sigma_f_lib(2);
sigs12   = sigma_s_lib(1, 1, 2);   % sigma_s (group 1 -> group 2)

%% ========================================================================
%  2. ANALYTICAL SOLUTION  (fundamental mode: m = n = 1)
% =========================================================================

m = 1;  n = 1;   % fundamental mode indices

% Geometric buckling (2D)
Bx2  = (m * pi / L)^2;
By2  = (n * pi / L)^2;
B2   = Bx2 + By2;

% Analytic eigenvalue lambda (= K_eff for fundamental mode)
k_analytic = (nusigf1 * (D2*B2 + siga2) + nusigf2 * sigs12) / ...
             ((D2*B2 + siga2) * (siga1 + sigs12 + D1*B2));

fprintf('==============================================\n');
fprintf('  2-GROUP 2D HOMOGENEOUS REACTOR — ANALYTICS  \n');
fprintf('==============================================\n\n');
fprintf('Geometric buckling B^2 = %.6e cm^-2\n', B2);
fprintf('Analytic K-eff         = %.10f\n\n',    k_analytic);

% Fast-to-thermal flux ratio (spatially constant)
ratio_1_to_2 = (D2*B2 + siga2) / sigs12;
fprintf('Flux ratio psi1/psi2   = %.6f\n\n', ratio_1_to_2);

% Normalisation constant k (2D adaptation of Eq. A.11)
%   integral of |sin(m*pi*x/L)*sin(n*pi*y/L)| over [0,L]x[0,L]
%   = (2L/(m*pi)) * (2L/(n*pi)) = 4*L^2 / (m*n*pi^2)
%   For m=n=1: 4*L^2/pi^2
integral_psi2 = 4 * L^2 / (m * n * pi^2);

% Normalisation condition (Eq. A.11 adapted to 2D, V_total = L^2):
%   1 = (1/L^2) * (sigf1*ratio_1_to_2 + sigf2) * k * integral_psi2
k_norm = L^2 / ((sigf1 * ratio_1_to_2 + sigf2) * integral_psi2);

fprintf('Normalisation constant k = %.6e\n\n', k_norm);

% Analytical flux function handles
phi2_func = @(x, y) k_norm ...
    .* sin(m * pi * x / L) ...
    .* sin(n * pi * y / L);

phi1_func = @(x, y) ratio_1_to_2 .* phi2_func(x, y);

%% ========================================================================
%  3. VERIFY ANALYTICAL FLUXES SATISFY THE DIFFUSION EQUATIONS
%     (numerical check via finite differences on a fine grid)
% =========================================================================

fprintf('--- Verification of analytical fluxes ---\n');
N_check = 500;
x_vec   = linspace(0, L, N_check+2);  x_vec = x_vec(2:end-1);
y_vec   = linspace(0, L, N_check+2);  y_vec = y_vec(2:end-1);
[X, Y]  = meshgrid(x_vec, y_vec);

P2 = phi2_func(X, Y);
P1 = phi1_func(X, Y);

% Laplacian of psi_2 analytically = -B^2 * psi_2
lap_P2_analytic = -B2 * P2;

% Residual of group-2 equation:
%   -D2*lap(psi2) + siga2*psi2 - sigs12*psi1 = 0
res2 = -D2 * lap_P2_analytic + siga2 * P2 - sigs12 * P1;
fprintf('Max |residual group 2| = %.2e  (should be ~0)\n', max(abs(res2(:))));

% Residual of group-1 equation:
%   -D1*lap(psi1) + (siga1+sigs12)*psi1 - (1/k)*nusigf2*psi2 - (1/k)*nusigf1*psi1 = 0
lap_P1_analytic = -B2 * P1;
res1 = -D1 * lap_P1_analytic + (siga1 + sigs12) * P1 ...
       - (1/k_analytic) * (nusigf1 * P1 + nusigf2 * P2);
fprintf('Max |residual group 1| = %.2e  (should be ~0)\n\n', max(abs(res1(:))));

%% ========================================================================
%  4. NUMERICAL CONVERGENCE STUDY
% =========================================================================

N_array       = [10, 20, 40, 80];
k_eff_num_vec = zeros(length(N_array), 1);

fprintf('%-10s  %-12s  %-14s  %-12s\n', 'Grid', 'Matrix size', 'k_numeric', 'err (pcm)');
fprintf('%s\n', repmat('-', 1, 54));

for i = 1:length(N_array)
    N = N_array(i);

    delta            = (L / N) * ones(N, 1);
    mesh             = Mesh_2D_FDM(delta, delta);
    region_materials = ones(N, N);   % (Ny x Nx) — homogeneous, all active

    nuclear_data = NuclearData_2D(region_materials, D_lib, sigma_a_lib, ...
                                  nu_sigf_lib, chi_lib, sigma_s_lib);

    solver = Solver_2D_FDM(mesh, nuclear_data);
    solver = solver.assembleMatrices();
    solver = solver.solveEigenvalues(1);

    k_eff_num_vec(i) = solver.keff(1);
    err_pcm          = abs(k_analytic - k_eff_num_vec(i)) * 1e5;

    fprintf('%-10s  %-12d  %-14.10f  %-12.2f\n', ...
        sprintf('%dx%d', N, N), N^2, k_eff_num_vec(i), err_pcm);
end

%% ========================================================================
%  5. PLOT ANALYTICAL FLUX SHAPES
% =========================================================================

x_plot = linspace(0, L, 200);
y_plot = linspace(0, L, 200);
[Xp, Yp] = meshgrid(x_plot, y_plot);

figure('Name', 'Analytical Fluxes — 2-Group 2D Homogeneous Reactor', 'Color', 'w');

subplot(1,2,1);
imagesc(x_plot, y_plot, phi1_func(Xp, Yp));
set(gca, 'YDir', 'normal'); colorbar; axis equal tight;
xlabel('X (cm)'); ylabel('Y (cm)');
title(sprintf('Fast flux \\psi_1(x,y)  [m=n=1]'));

subplot(1,2,2);
imagesc(x_plot, y_plot, phi2_func(Xp, Yp));
set(gca, 'YDir', 'normal'); colorbar; axis equal tight;
xlabel('X (cm)'); ylabel('Y (cm)');
title(sprintf('Thermal flux \\psi_2(x,y)  [m=n=1]'));

sgtitle(sprintf('K_{eff} = %.5f', k_analytic));

%% ========================================================================
%  6. SAVE RESULTS
% =========================================================================

save_path = fullfile(fileparts(mfilename('fullpath')), '..', 'data', ...
                     'homo_2EG_2D_exact_results.mat');

save(save_path, 'N_array', 'k_eff_num_vec', 'k_analytic', ...
     'phi1_func', 'phi2_func', 'k_norm', 'ratio_1_to_2', 'B2', ...
     'D_lib', 'sigma_a_lib', 'nu_sigf_lib', 'chi_lib', 'sigma_s_lib', 'L');

fprintf('\nResults saved to: %s\n', save_path);