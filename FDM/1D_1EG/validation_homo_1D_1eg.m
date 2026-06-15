% =========================================================================
% SCRIPT: Comparison FDM (Ruben) vs FDM (Jaime) vs FEM Degree 3 (Jaime)
% =========================================================================
addpath('../..') 
clearvars; close all; clc;

%% 1. PROBLEM PARAMETERS (1G Homogeneous)
L = 350;
FEM_DEGREE = 3;

% Array of different mesh sizes (Degrees of Freedom) to test
DOF_array = [43, 100, 500, 1000]; 

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

% START LOOP OVER MESH SIZES
for target_dof = DOF_array

    fprintf('--- STARTING VALIDATION FOR DOF = %d ---\n\n', target_dof);

    %% 2. SOLUTION A: MESH-CENTERED FDM (RUBEN)
    N_Ruben = target_dof; % cells
    meshR = Mesh_1D_FDM(L, N_Ruben);
    matsR = Materials(1, D_val, SigA_val, NuSigF_val);
    probR = Solver_1D_1EG_FDM(meshR, matsR);
    probR = probR.assembleMatrices().solveEigenvalues(1);
    k_eff_R = probR.keff(1);
    size_R = size(probR.A,1);
    % Extract flux (including boundaries for the plot)
    phi_num_R = probR.phi(:, 1);
    x_R = probR.full_mesh_coordinates;
    % Normalize the analytic flux (numeric flux is already normalized)
    phi_ana_R = phi_analytic_func(x_R);
    phi_ana_R = phi_ana_R / max(abs(phi_ana_R));
    % Calculate Absolute RMSE
    rmse_R = rmse(phi_ana_R, phi_num_R);

    %% 3. SOLUTION B: POINT-CENTERED FDM (JAIME)
    N_J_FDM = target_dof; % points
    cell_sizes_fdm = (L/N_J_FDM) * ones(1, N_J_FDM);
    cell_sizes_fdm(end) = L - sum(cell_sizes_fdm(1:end-1));
    material_map_fdm = 2 * ones(1, N_J_FDM);
    % Initialize Jaime's classes for FDM
    meshJ_FDM = Malla1D(L, N_J_FDM, 1, cell_sizes_fdm);
    materialsJ_FDM = Materiales1D1g(material_map_fdm, D_array, SigA_array, NuSigF_array, meshJ_FDM);
    probJ_FDM = ProblemaDifusion1D1g(meshJ_FDM, materialsJ_FDM, []);
    probJ_FDM = probJ_FDM.ensamblar_matrices_fdm().aplicar_cc().resolver_autovalor(1);
    k_eff_J_FDM = probJ_FDM.keff(1);
    size_J_FDM = size(probJ_FDM.A,1);
    % Extract flux and generate matching point-centered coordinates
    phi_num_J_FDM = probJ_FDM.phi_inc(2:end-1, 1);
    x_full_fdm = linspace(0, L, N_J_FDM)'; 
    x_J_FDM = x_full_fdm(2:end-1); 
    % Error Logic
    phi_ana_J_FDM = phi_analytic_func(x_J_FDM);
    phi_ana_J_FDM = phi_ana_J_FDM / max(abs(phi_ana_J_FDM));
    rmse_J_FDM = rmse(phi_ana_J_FDM, phi_num_J_FDM);

    %% 4. SOLUTION C: FEM DEGREE 3 (JAIME)
    % Element calculation: dynamically computed based on target_dof
    N_FEM = max(1, round((target_dof - 1) / FEM_DEGREE)); 
    cell_sizes = (L/N_FEM) * ones(1, N_FEM);
    cell_sizes(end) = L - sum(cell_sizes(1:end-1)); % This is to prevent a size error because of the N_FEM approximation
    material_map = 2 * ones(1, N_FEM); % Set all to Material 2
    % Initialize Jaime's classes for FEM
    meshJ_FEM = Malla1D(L, N_FEM, FEM_DEGREE, cell_sizes);
    materialsJ_FEM = Materiales1D1g(material_map, D_array, SigA_array, NuSigF_array, meshJ_FEM);
    elementJ = ElementoFinito(FEM_DEGREE);
    probJ_FEM = ProblemaDifusion1D1g(meshJ_FEM, materialsJ_FEM, elementJ);
    probJ_FEM = probJ_FEM.ensamblar_matrices_fem().aplicar_cc().resolver_autovalor(1);
    k_eff_F = probJ_FEM.keff(1);
    size_F = size(probJ_FEM.A, 1);
    % Extract flux (Jaime uses the full vector, including boundaries)
    phi_num_F = probJ_FEM.phi_inc(:, 1);
    x_F = meshJ_FEM.x_nodos(:);
    % Error Logic
    phi_ana_F = phi_analytic_func(x_F);
    phi_num_F = phi_num_F / max(abs(phi_num_F));
    %phi_ana_F = phi_ana_F / max(abs(phi_ana_F));
    rmse_F = rmse(phi_ana_F, phi_num_F);

    %% 5. PRINT RESULTS
    err_pcm_R = abs(k_analytic - k_eff_R) * 1e5;
    err_pcm_J_FDM = abs(k_analytic - k_eff_J_FDM) * 1e5;
    err_pcm_F = abs(k_analytic - k_eff_F) * 1e5;
    fprintf('========================================================================================\n');
    fprintf('                           DIRECT COMPARISON (DOF = %d)                                 \n', target_dof);
    fprintf('========================================================================================\n');
    fprintf('Analytic K-eff: %.10f\n', k_analytic);
    fprintf('----------------------------------------------------------------------------------------\n');
    fprintf('%-18s | %-18s | %-18s | %-20s\n', 'Metric', 'FDM (Ruben)', 'FDM (Jaime)', 'FEM Deg 3 (Jaime)');
    fprintf('-------------------|--------------------|--------------------|------------------------\n');
    fprintf('%-18s | %-18d | %-18d | %-20d\n', 'Matrix Size', size_R, size_J_FDM, size_F);
    fprintf('%-18s | %-18.8f | %-18.8f | %-20.8f\n', 'Numeric K-eff', k_eff_R, k_eff_J_FDM, k_eff_F);
    fprintf('%-18s | %-18.3f | %-18.3f | %-20.3f\n', 'K-eff Err (pcm)', err_pcm_R, err_pcm_J_FDM, err_pcm_F);
    fprintf('%-18s | %-18.6e | %-18.6e | %-20.6e\n', 'Absolute RMSE', rmse_R, rmse_J_FDM, rmse_F);
    fprintf('========================================================================================\n');

    %% 6. COMPARATIVE PLOTS
    figure('Name', sprintf('Validation All Methods DOF=%d', target_dof), 'Color', 'w', 'Position', [100, 100, 800, 700]);
    % --- Normalized Flux ---
    hold on; grid on;
    plot(x_R, phi_ana_R, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytic Solution');
    plot(x_R, phi_num_R, 'b--o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'FDM (Ruben)');
    plot(x_J_FDM, phi_num_J_FDM, 'g--s', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'FDM (Jaime)');
    plot(x_F, phi_num_F, 'r--x', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'FEM Deg 3 (Jaime)');
    xlabel('Position x (cm)');
    ylabel('Normalized Flux \phi(x)');
    title(sprintf('Flux Shape Comparison (1G Homogeneous, DOF=%d)', target_dof));
    legend('Location', 'best');

end