%% homo_1EG_example
% This script extracts the validated solution for a Homogeneous 1EG problem
% and saves it in a file called "homo_1EG_exact_results" in order
% to define a test script

clearvars; close all; clc;
addpath("../")

% PROBLEM PARAMETERS (1G Homogeneous)
L = 350;

% Array of different mesh sizes (Degrees of Freedom or matrix sizes)
N_array = [43, 100, 500, 1000]; 

% Physical properties (Index 2 = Homogeneous Fuel)
D_array       = [1.446, 0.776];
SigA_array    = [0.0077, 0.0244];
NuSigF_array  = [0.0000, 0.0260];
D_val = D_array(2);
SigA_val = SigA_array(2);
NuSigF_val = NuSigF_array(2);

% Exact Analytical Solution (Mode 1)
k_analytic = NuSigF_val / (D_val * (pi / L)^2 + SigA_val);
phi_analytic_func = @(x) sin(pi * x / L);

disp("Processing the problem...")
fprintf("Analytic K-eff: %.10f\n", k_analytic)
disp("Numerical values:")

% Create vector for saving the numerical k_eff values
k_eff_num_vec = zeros(length(N_array));

% Compute numerical solution for each mesh size
for i = 1:length(N_array)
    N = N_array(i);
    mesh = Mesh_1D_FDM(L, N);
    materials = Materials(1, D_val, SigA_val, NuSigF_val);
    solver = Solver_1D_1EG_FDM(mesh, materials);
    solver = solver.assembleMatrices().solveEigenvalues(1);
    
    % Extract the computed k_eff
    k_eff_num = solver.keff(1);
    
    % Save the computed k_eff and print on screen
    k_eff_num_vec(i) = solver.keff(1);
    fprintf("Matrix size: %d | Numeric k_eff: %.10f\n", N, k_eff_num)
end

% Save the variable values in a file
save("homo_1EG_exact_results.mat", "N_array", "k_eff_num_vec", "k_analytic");
