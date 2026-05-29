%% het_1EG_example
% This script extracts the validated solution for a Heterogeneous
% 1EG problem and saves it in a file called "het_1EG_exact_results"
% in order to define a test script

clearvars; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')); 

% PROBLEM PARAMETERS (1G Heterogeneous)
L = 350;
x_int1 = 25;   % Coordinate of Interface 1 (Reflector -> Core)
x_int2 = 325;  % Coordinate of Interface 2 (Core -> Reflector)

% Array of different mesh sizes (Degrees of Freedom) to test
N_array = [43, 100, 500, 1000];

% Physical properties
% Index 1 = Reflector, Index 2 = Homogeneous Fuel
D_lib       = [1.446; 0.776];
sigma_a_lib = [0.0077; 0.0244];
nu_sigf_lib = [0.0000; 0.0260];
chi_lib     = [0.0; 1.0]; 
sigma_s_lib = zeros(2, 1, 1);

% EXACT ANALYTICAL SOLUTION (Mode 1)

% Map properties
D_1 = D_lib(1);
Sig_a1 = sigma_a_lib(1);
D_2 = D_lib(2);
Sig_a2 = sigma_a_lib(2);
nuSig_f = nu_sigf_lib(2);

% Explicitly define the physical dimensions for the symmetric solver
a_core = (x_int2 - x_int1) / 2; % core's half length (150 cm)
b_refl = x_int1;                % reflector's length (25 cm)
x_center = x_int1 + a_core;     % center of the reactor (175 cm)
kappa_1 = sqrt(Sig_a1 / D_1);

% Calculate Initial Guess (Bare Reactor)
d_ext = 2.13 * D_2;
a_ext = a_core + d_ext;
Bg = pi / (2.0 * a_ext);
k_bare = nuSig_f / (D_2 * Bg^2 + Sig_a2);

% Calculate infinite multiplication factor (upper physical limit)
k_inf = nuSig_f / Sig_a2;

% Solve Criticality Condition for k_eff using physical bracketing
options = optimset('Display', 'off', 'TolX', 1e-12);
k_analytic = fzero(@(k) criticality_eq(k, nuSig_f, Sig_a2, D_2, D_1, kappa_1, a_core, b_refl), [k_bare, k_inf - 1e-10], options);

% Calculate Final Buckling and Flux Constants
B_2 = sqrt((nuSig_f / k_analytic - Sig_a2) / D_2);
A_const = 1.0; % Max flux at center (arbitrary)
C_const = A_const * cos(B_2 * a_core) / sinh(kappa_1 * b_refl);

% Function to evaluate exact analytical flux using the symmetric formulation
phi_analytic_func = @(x) eval_analytical_flux_symmetric(x, x_center, a_core, b_refl, B_2, kappa_1, A_const, C_const);

disp("Processing the problem...")
fprintf('Analytic K-eff: %.10f\n', k_analytic);

disp("Numerical values:")

size_vec = zeros(length(N_array));

phi_analytic_array = {};

% Compute numerical solution for each mesh size
for i = 1:length(N_array)  
    N = N_array(i);
    
    region_lengths = [b_refl; 2 * a_core; L - x_int2];
    region_materials = [1; 2; 1];
    cells_per_region = ceil(N / 3) * ones(3,1);
    
    mesh = Mesh_1D_FDM(region_lengths, cells_per_region);

    full_mesh_coordinates = [0; mesh.cell_centers; 0];
        
    k_eff_num_vec(i) = k_eff_num;
    
    phi_analytic_array{i} = phi_analytic_func(full_mesh_coordinates) / max(phi_analytic_func(full_mesh_coordinates));

    size_vec(i) = size(solver.A,1);
    fprintf("Matrix size: %d | Numeric k_eff: %.10f\n", size_vec(i), k_eff_num)
end

% Save the variable values in a file
save("../data/het_1EG_exact_results.mat", "size_vec", "k_analytic", "phi_analytic_array");

% =========================================================================
% LOCAL FUNCTIONS FOR ANALYTICAL SOLUTION
% =========================================================================
function res = criticality_eq(k, nuSig_f, Sig_a2, D_2, D_1, kappa_1, a_core, b_refl)
    B_2_sq = (nuSig_f / k - Sig_a2) / D_2;
    if B_2_sq <= 0
        res = 1e6; 
    else
        B_2 = sqrt(B_2_sq);
        coth_kb = 1.0 / tanh(kappa_1 * b_refl);
        res = D_2 * B_2 * tan(B_2 * a_core) - D_1 * kappa_1 * coth_kb;
    end
end

function phi_val = eval_analytical_flux_symmetric(x_array, x_center, a_core, b_refl, B_2, kappa_1, A_const, C_const)
    phi_val = zeros(size(x_array));
    for i = 1:length(x_array)
        x_shifted = abs(x_array(i) - x_center);
        if x_shifted <= a_core
            phi_val(i) = A_const * cos(B_2 * x_shifted);
        else
            phi_val(i) = C_const * sinh(kappa_1 * (a_core + b_refl - x_shifted));
        end
    end
end