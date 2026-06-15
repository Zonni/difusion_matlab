% =========================================================================
% SCRIPT: Comparison FDM (Ruben) vs FDM (Jaime) vs FEM Degree 3 (Jaime)
% PROBLEM: 1G Heterogeneous (Reflector - Core - Reflector)
% =========================================================================
addpath('../..') 
clearvars; close all; clc;

%% 1. PROBLEM PARAMETERS (1G Heterogeneous)
L = 350;
x_int1 = 25;   % Coordinate of Interface 1 (Reflector -> Core)
x_int2 = 325;  % Coordinate of Interface 2 (Core -> Reflector)
FEM_DEGREE = 4;

% Array of different mesh sizes (Degrees of Freedom) to test
DOF_array = [43, 100, 500, 1000];

% Physical properties
% Index 1 = Reflector, Index 2 = Homogeneous Fuel
D_array       = [1.446, 0.776];
SigA_array    = [0.0077, 0.0244];
NuSigF_array  = [0.0000, 0.0260];

%% 2. EXACT ANALYTICAL SOLUTION
fprintf('--- CALCULATING EXACT ANALYTICAL SOLUTION ---\n');

% Map properties
D_1 = D_array(1);
Sig_a1 = SigA_array(1);
D_2 = D_array(2);
Sig_a2 = SigA_array(2);
nuSig_f = NuSigF_array(2);

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

fprintf('Initial guess (Bare k_eff) : %.5f\n', k_bare);
fprintf('Analytic K-eff: %.10f\n\n', k_analytic);

rmse_R_vec = zeros(length(DOF_array),1);
rmse_F_vec = zeros(length(DOF_array),1);
i = 1;

%% START LOOP OVER MESH SIZES
for target_dof = DOF_array
    fprintf('--- STARTING VALIDATION FOR DOF = %d ---\n', target_dof);
    
    %% 3. SOLUTION A: MESH-CENTERED FDM (RUBEN)
    N_Ruben = target_dof; % cells
    
    region_lengths = [b_refl; 2 * a_core; L - x_int2];
    region_materials = [1; 2; 1];
    cells_per_region = ceil(target_dof / 3) * ones(3,1);
    
    meshR = Mesh_1D_FDM(region_lengths, cells_per_region);
    
    matsR = Materials(region_materials, D_array, SigA_array, NuSigF_array); 
    probR = Solver_1D_1EG_FDM(meshR, matsR);
    probR = probR.assembleMatrices().solveEigenvalues(1);
    
    k_eff_R = probR.keff(1);
    size_R = size(probR.A,1);
    
    % Extract flux
    phi_num_R = probR.phi(:, 1);
    x_R = probR.full_mesh_coordinates;
    
    % Normalize analytic
    phi_ana_R = phi_analytic_func(x_R);
    phi_ana_R = phi_ana_R / max(abs(phi_ana_R)); 
    
    rmse_R = rmse(phi_ana_R, phi_num_R);
    rmse_R_vec(i) = rmse_R;
    
    %% 4. SOLUTION B: POINT-CENTERED FDM (JAIME)
    N_J_FDM = target_dof; % points
    cell_sizes_fdm = (L/N_J_FDM) * ones(1, N_J_FDM);
    cell_sizes_fdm(end) = L - sum(cell_sizes_fdm(1:end-1)); 
    
    % Create material map based on cell centers
    x_centers_J = cumsum(cell_sizes_fdm) - cell_sizes_fdm/2;
    material_map_fdm = ones(1, N_J_FDM);
    material_map_fdm(x_centers_J >= x_int1 & x_centers_J <= x_int2) = 2;
    
    % Initialize Jaime's classes for FDM
    meshJ_FDM = Malla1D(L, N_J_FDM, 1, cell_sizes_fdm);
    materialsJ_FDM = Materiales1D1g(material_map_fdm, D_array, SigA_array, NuSigF_array, meshJ_FDM);
    probJ_FDM = ProblemaDifusion1D1g(meshJ_FDM, materialsJ_FDM, []);
    probJ_FDM = probJ_FDM.ensamblar_matrices_fdm().aplicar_cc().resolver_autovalor(1);
    
    k_eff_J_FDM = probJ_FDM.keff(1);
    size_J_FDM = size(probJ_FDM.A,1);
    
    % Extract flux
    phi_num_J_FDM = probJ_FDM.phi_inc(2:end-1, 1);
    x_full_fdm = linspace(0, L, N_J_FDM)'; 
    x_J_FDM = x_full_fdm(2:end-1); 
    
    % Normalize analytic
    phi_ana_J_FDM = phi_analytic_func(x_J_FDM);
    phi_ana_J_FDM = phi_ana_J_FDM / max(abs(phi_ana_J_FDM));
    
    rmse_J_FDM = rmse(phi_ana_J_FDM, phi_num_J_FDM);
    
    %% 5. SOLUTION C: FEM DEGREE 4 (JAIME)
    N_FEM = max(1, round((target_dof - 1) / FEM_DEGREE)); 
    cell_sizes = (L/N_FEM) * ones(1, N_FEM);
    cell_sizes(end) = L - sum(cell_sizes(1:end-1)); 
    
    % Create material map based on element centers
    x_centers_FEM = cumsum(cell_sizes) - cell_sizes/2;
    material_map_FEM = ones(1, N_FEM);
    material_map_FEM(x_centers_FEM >= x_int1 & x_centers_FEM <= x_int2) = 2;
    
    % Initialize Jaime's classes for FEM
    meshJ_FEM = Malla1D(L, N_FEM, FEM_DEGREE, cell_sizes);
    materialsJ_FEM = Materiales1D1g(material_map_FEM, D_array, SigA_array, NuSigF_array, meshJ_FEM);
    elementJ = ElementoFinito(FEM_DEGREE);
    probJ_FEM = ProblemaDifusion1D1g(meshJ_FEM, materialsJ_FEM, elementJ);
    probJ_FEM = probJ_FEM.ensamblar_matrices_fem().aplicar_cc().resolver_autovalor(1);
    
    k_eff_F = probJ_FEM.keff(1);
    size_F = size(probJ_FEM.A, 1);
    
    % Extract flux 
    phi_num_F = probJ_FEM.phi_inc(:, 1);
    x_F = meshJ_FEM.x_nodos(:);
    
    % Normalize analytic
    phi_ana_F = phi_analytic_func(x_F);
    phi_ana_F = phi_ana_F / max(abs(phi_ana_F));
    
    rmse_F = rmse(phi_ana_F, phi_num_F);
    rmse_F_vec(i) = rmse_F;
    
    %% 6. PRINT RESULTS
    err_pcm_R = abs(k_analytic - k_eff_R) * 1e5;
    err_pcm_J_FDM = abs(k_analytic - k_eff_J_FDM) * 1e5;
    err_pcm_F = abs(k_analytic - k_eff_F) * 1e5;
    
    fprintf('========================================================================================\n');
    fprintf('                      DIRECT COMPARISON (DOF = %d) - HETEROGENEOUS                      \n', target_dof);
    fprintf('========================================================================================\n');
    fprintf('%-18s | %-18s | %-18s | %-20s\n', 'Metric', 'FDM (Ruben)', 'FDM (Jaime)', 'FEM Deg 4 (Jaime)');
    fprintf('-------------------|--------------------|--------------------|------------------------\n');
    fprintf('%-18s | %-18d | %-18d | %-20d\n', 'Matrix Size', size_R, size_J_FDM, size_F);
    fprintf('%-18s | %-18.8f | %-18.8f | %-20.8f\n', 'Numeric K-eff', k_eff_R, k_eff_J_FDM, k_eff_F);
    fprintf('%-18s | %-18.3f | %-18.3f | %-20.3f\n', 'K-eff Err (pcm)', err_pcm_R, err_pcm_J_FDM, err_pcm_F);
    fprintf('%-18s | %-18.6e | %-18.6e | %-20.6e\n', 'Absolute RMSE', rmse_R, rmse_J_FDM, rmse_F);
    fprintf('========================================================================================\n\n');
    
    %% 7. COMPARATIVE PLOTS
    figure('Name', sprintf('Validation Heterogeneous DOF=%d', target_dof), 'Color', 'w', 'Position', [100, 100, 800, 700]);
    hold on; grid on;
    
    x_ana = linspace(0, 350, 100);
    phi_ana = phi_analytic_func(x_ana);
    phi_ana = phi_ana / max(phi_ana);

    plot(x_ana, phi_ana, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytic Solution');
    plot(x_R, phi_num_R, 'b--o', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'FDM (Ruben)');
    plot(x_J_FDM, phi_num_J_FDM, 'g--s', 'LineWidth', 1.5, 'MarkerSize', 5, 'DisplayName', 'FDM (Jaime)');
    plot(x_F, phi_num_F, 'r--x', 'LineWidth', 1.5, 'MarkerSize', 6, 'DisplayName', 'FEM Deg 4 (Jaime)');
    
    xline(x_int1, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xline(x_int2, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    xlabel('Position x (cm)');
    ylabel('Normalized Flux \phi(x)');
    title(sprintf('Flux Shape Comparison (1G Heterogeneous, DOF=%d)', target_dof));
    legend('Location', 'best');

    i = i + 1;
end

figure('Color', 'w', 'Name', 'Convergence Study');
loglog(DOF_array, rmse_R_vec, '-ok', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', 'Ruben FDM');
grid on; hold on;
loglog(DOF_array, rmse_F_vec, '-sr', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'Jaime FEM (Deg 4)'); % NEW: FEM Plot

% Add labels and title
xlabel('Number of Nodes (DOF)');
ylabel('Flux RMSE (Absolute)');
title('Convergence Comparison: FDM vs FEM');

ref_line = (rmse_R_vec(1)) * (DOF_array(1)./DOF_array).^2;
loglog(DOF_array, ref_line, '--k', 'DisplayName', 'Theoretical O(h^2)');
legend('Location', 'southwest');

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