%% nakata_martin_benchmark.m 
% Solves the 2D two-group benchmark from: 
%   Nakata, H. and Martin, W. R. 
%   "The Finite Element Response Matrix Method" 
%   Nuclear Science and Engineering: 85, 289-305 (1983) 
% 
% The reactor is an irregular octagonal shape on a 17x17 assembly grid. 
% Each assembly is 23.1226 cm x 23.1226 cm. 
% 
% A mesh refinement study is performed by subdividing each assembly into 
% N x N uniform sub-cells (N = 1, 2, 3, 4, 5, 6).  For each level the script 
% records DOF, assembly time, solve time, k-eff, and error in pcm. 

clearvars; close all; clc; 
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', '2D')); 

%% ======================================================================== 
%  1. NUCLEAR DATA  (fixed across all refinement levels) 
% ========================================================================= 
D_lib = [ ... 
    1.4360,  0.3635 ; 
    1.4366,  0.3636 ; 
    1.3200,  0.2772 ; 
    1.4389,  0.3638 ; 
    1.4381,  0.3665 ; 
    1.4385,  0.3665 ; 
    1.4389,  0.3679 ; 
    1.4393,  0.3680 ];

sigma_a_lib = [ ... 
    0.0095042,  0.0750580 ; 
    0.0096785,  0.0784360 ; 
    0.0026562,  0.0715960 ; 
    0.0103630,  0.0914080 ; 
    0.0100030,  0.0848280 ; 
    0.0101320,  0.0873140 ; 
    0.0101650,  0.0880240 ; 
    0.0102940,  0.0905100 ];

nu_sigf_lib = [ ... 
    0.0058708,  0.0960670 ; 
    0.0061908,  0.1035800 ; 
    0.0000000,  0.0000000 ; 
    0.0074527,  0.1323600 ; 
    0.0061908,  0.1035800 ; 
    0.0064285,  0.1091100 ; 
    0.0061908,  0.1035800 ; 
    0.0064285,  0.1091100 ];

chi_lib = repmat([1.0, 0.0], 8, 1);

Sigma_12             = [0.017754; 0.017621; 0.023106; 0.017101; ... 
                        0.017290; 0.017192; 0.017125; 0.017027];

sigma_s_lib          = zeros(8, 2, 2); 

sigma_s_lib(:, 1, 2) = Sigma_12; 

%% ======================================================================== 
%  2. ASSEMBLY-LEVEL MATERIAL GRID  (17 x 17) 
%  0 = outside reactor boundary 
% ========================================================================= 
assembly_mat = zeros(17, 17); 
assembly_mat( 1,:) = [0 0 0 0 3 3 3 3 3 3 3 3 3 0 0 0 0]; 
assembly_mat( 2,:) = [0 0 3 3 3 4 4 4 4 4 4 4 3 3 3 0 0]; 
assembly_mat( 3,:) = [0 3 3 4 4 8 1 1 1 1 1 8 4 4 3 3 0]; 
assembly_mat( 4,:) = [0 3 4 4 5 1 7 1 7 1 7 1 5 4 4 3 0]; 
assembly_mat( 5,:) = [3 3 4 5 2 8 2 8 1 8 2 8 2 5 4 3 3]; 
assembly_mat( 6,:) = [3 4 8 1 8 2 8 2 6 2 8 2 8 1 8 4 3]; 
assembly_mat( 7,:) = [3 4 1 7 2 8 1 8 2 8 1 8 2 7 1 4 3]; 
assembly_mat( 8,:) = [3 4 1 1 8 2 8 1 8 1 8 2 8 1 1 4 3]; 
assembly_mat( 9,:) = [3 4 1 7 1 6 2 8 1 8 2 6 1 7 1 4 3]; 
assembly_mat(10,:) = [3 4 1 1 8 2 8 1 8 1 8 2 8 1 1 4 3]; 
assembly_mat(11,:) = [3 4 1 7 2 8 1 8 2 8 1 8 2 7 1 4 3]; 
assembly_mat(12,:) = [3 4 8 1 8 2 8 2 6 2 8 2 8 1 8 4 3]; 
assembly_mat(13,:) = [3 3 4 5 2 8 2 8 1 8 2 8 2 5 4 3 3]; 
assembly_mat(14,:) = [0 3 4 4 5 1 7 1 7 1 7 1 5 4 4 3 0]; 
assembly_mat(15,:) = [0 3 3 4 4 8 1 1 1 1 1 8 4 4 3 3 0]; 
assembly_mat(16,:) = [0 0 3 3 3 4 4 4 4 4 4 4 3 3 3 0 0]; 
assembly_mat(17,:) = [0 0 0 0 3 3 3 3 3 3 3 3 3 0 0 0 0]; 

assembly_size = 23.1226;   % cm 

%% ======================================================================== 
%  3. REFERENCE K-EFF 
%  Nakata & Martin (1983), Table II. 
% ========================================================================= 
k_ref = 1.025110; 

%% ======================================================================== 
%  4. REFINEMENT STUDY 
% ========================================================================= 
N_refine  = [1, 2, 3, 4, 5, 6]; % sub-cells per assembly side 
n_levels  = length(N_refine); 

keff_vec = zeros(n_levels, 1); 
dof_vec = zeros(n_levels, 1); 
t_asm_vec = zeros(n_levels, 1); 
t_sol_vec = zeros(n_levels, 1); 
err_pcm = zeros(n_levels, 1); 
 
power_history_17x17 = zeros(17, 17, 2);  % Store coarse and fine normalized maps

fprintf('Running refinement study...\n\n'); 

for refinment = 1:n_levels 
    N = N_refine(refinment); 
     
    % Expand assembly grid 
    region_materials = kron(assembly_mat, ones(N, N)); 
     
    % Mesh 
    region_lengths   = assembly_size * ones(17, 1); 
    cells_per_region = N * ones(17, 1); 
    mesh = Mesh_2D_FDM(region_lengths, cells_per_region, ... 
                       region_lengths, cells_per_region); 
                    
    nuclear_data = NuclearData_2D(region_materials, D_lib, sigma_a_lib, ... 
                                  nu_sigf_lib, chi_lib, sigma_s_lib); 
                               
    % Assembly timing 
    solver = Solver_2D_FDM(mesh, nuclear_data, 'vacuum'); 
    t0 = tic; 
    solver = solver.assembleMatrices(); 
    t_asm_vec(refinment) = toc(t0); 
     
    % Solve timing 
    t0 = tic; 
    solver = solver.solveEigenvalues(1); 
    t_sol_vec(refinment) = toc(t0); 
     
    keff_vec(refinment) = solver.keff(1); 
    dof_vec(refinment) = 2 * numel(solver.active_cells);    
    err_pcm(refinment) = abs(keff_vec(refinment) - k_ref) * 1e5; 
     
    fprintf('N=%d | grid %dx%d | DOF=%d | k=%.6f | err=%.1f pcm | t=%.2fs\n', ... 
            N, 17*N, 17*N, dof_vec(refinment), keff_vec(refinment), err_pcm(refinment), ... 
            t_asm_vec(refinment) + t_sol_vec(refinment)); 
     
    % --------------------------------------------------------------------- 
    % POWER COMPUTATION & COLLAPSE (Only for first and last refinement) 
    % --------------------------------------------------------------------- 
    if refinment == 1 || refinment == n_levels 
        
        % 1. Get the raw fine-mesh power
        solver = solver.computePower(1); 
        raw_power_density = solver.power(:, :, 1); 
         
        % 2. Collapse the raw power to the 17x17 assembly level 
        level_power = zeros(17, 17); 
        for row = 1:17 
            for col = 1:17 
                row_idx = (row - 1) * N + 1 : row * N; 
                col_idx = (col - 1) * N + 1 : col * N; 
                block_power = raw_power_density(row_idx, col_idx); 
                level_power(row, col) = mean(block_power(:), 'omitnan'); 
            end 
        end 
        
        % 3. Normalize at the ASSEMBLY level (Matches your validation report exactly)
        active_assemblies = (assembly_mat > 0); 
        avg_core_power = mean(level_power(active_assemblies)); 
        level_power_norm = level_power ./ avg_core_power; 
         
        % Store in the preallocated history tensor 
        if refinment == 1 
            power_history_17x17(:, :, 1) = level_power_norm; 
        else 
            power_history_17x17(:, :, 2) = level_power_norm; 
        end 
    end 
    
    last_solver = solver;
end 

%% ======================================================================== 
%  5. RESULTS TABLE 
% ========================================================================= 
fprintf('\n'); 
sep = repmat('=', 1, 94); 
fprintf('%s\n', sep); 
fprintf('  Nakata-Martin (1983) — Mesh Refinement Study\n'); 
fprintf('  Reference k-eff: %.6f\n', k_ref); 
fprintf('%s\n', sep); 
fprintf('%-6s  %-10s  %-8s  %-12s  %-12s  %-10s  %-10s  %-10s\n', ... 
        'N', 'Grid', 'DOF', 'k_eff', 'err (pcm)', 'Order (p)', 't_asm (s)', 't_sol (s)'); 
fprintf('%s\n', repmat('-', 1, 94)); 
for refinment = 1:n_levels 
    N = N_refine(refinment); 
     
    % Calculate empirical order of convergence (p) 
    if refinment == 1 
        order_str = '-'; 
    else 
        % Since h is proportional to 1/N, the ratio h_{i-1}/h_i is N_i / N_{i-1} 
        p = log(err_pcm(refinment-1) / err_pcm(refinment)) / log(N_refine(refinment) / N_refine(refinment-1)); 
        order_str = sprintf('%.2f', p); 
    end 
     
    fprintf('%-6d  %-10s  %-8d  %-12.6f  %-12.1f  %-10s  %-10.3f  %-10.3f\n', ... 
            N, sprintf('%dx%d', 17*N, 17*N), dof_vec(refinment), ... 
            keff_vec(refinment), err_pcm(refinment), order_str, t_asm_vec(refinment), t_sol_vec(refinment)); 
end 
fprintf('%s\n', sep); 

%% ======================================================================== 
%  6. CONVERGENCE PLOTS 
% ========================================================================= 
figure('Name', 'Nakata-Martin Convergence Study', 'Color', 'w', ... 
       'Position', [80, 80, 1200, 400]); 
% --- Panel 1: k-eff error vs DOF (log-log) with O(h^2) reference -------- 
subplot(1, 3, 1); 
h2_ref = err_pcm(1) * (dof_vec(1) ./ dof_vec); 
loglog(dof_vec, err_pcm, 'o-b', 'LineWidth', 1.8, 'MarkerSize', 7, ... 
       'MarkerFaceColor', 'b'); 
hold on; 
loglog(dof_vec, h2_ref, '--k', 'LineWidth', 1.2); 
grid on; box on; 
legend('Numerical', 'O(h^2) reference', 'Location', 'southwest', ... 
       'FontSize', 9); 
xlabel('Degrees of Freedom'); 
ylabel('|k_{eff} - k_{ref}|  (pcm)'); 
title('k-eff Convergence'); 

% --- Panel 2: CPU time breakdown vs DOF (log-log) ------------------------ 
subplot(1, 3, 2); 
loglog(dof_vec, t_asm_vec, 's--r', 'LineWidth', 1.8, 'MarkerSize', 7, ... 
       'MarkerFaceColor', 'r'); 
hold on; 
loglog(dof_vec, t_sol_vec, '^--g', 'LineWidth', 1.8, 'MarkerSize', 7, ... 
       'MarkerFaceColor', 'g'); 
loglog(dof_vec, t_asm_vec + t_sol_vec, 'o-k', 'LineWidth', 1.8, ... 
       'MarkerSize', 7, 'MarkerFaceColor', 'k'); 
grid on; box on; 
legend('Assembly', 'Solve', 'Total', 'Location', 'northwest', 'FontSize', 9); 
xlabel('Degrees of Freedom'); 
ylabel('CPU time (s)'); 
title('Computational Cost'); 

% --- Panel 3: k-eff value vs refinement level ---------------------------- 
subplot(1, 3, 3); 
plot(N_refine, keff_vec, 'o-b', 'LineWidth', 1.8, 'MarkerSize', 7, ... 
     'MarkerFaceColor', 'b'); 
hold on; 
yline(k_ref, '--r', 'LineWidth', 1.5); 
grid on; box on; 
legend('Numerical k_{eff}', sprintf('Reference  k = %.6f', k_ref), ... 
       'Location', 'best', 'FontSize', 9); 
xlabel('Sub-cells per assembly side  N'); 
ylabel('k_{eff}'); 
xticks(N_refine); 
title('k-eff vs Mesh Level'); 
sgtitle('Nakata-Martin (1983) — Two-Group 2D Benchmark', 'FontSize', 12); 

%% ======================================================================== 
%  7. FLUX PLOT AT FINEST LEVEL 
% ========================================================================= 
last_solver.plotPhi(1); 

%% ======================================================================== 
%  8. PLOT ASSEMBLY POWER DISTRIBUTION 
% ========================================================================= 
fprintf('Generating relative assembly power plot for the finest mesh...\n'); 

% Extract the normalized power we already calculated in the loop
relative_power = power_history_17x17(:, :, 2); 
active_assemblies = (assembly_mat > 0); 

% Mask the reflector/outside nodes with NaN so they appear blank on the plot 
relative_power(~active_assemblies) = NaN; 

% Create the graphic 
figure('Name', 'BIBLIS Relative Assembly Power', 'Color', 'w', ... 
       'Position', [100, 100, 700, 600]); 
    
% Plot a 2D heatmap  
imagesc(relative_power, 'AlphaData', ~isnan(relative_power)); 
axis image;          % Force square pixels 
colormap jet;    % Standard MATLAB colormap 
set(gca, 'Color', [0.9 0.9 0.9]); % Light grey background for empty regions 

% Add colorbar and labels 
c = colorbar; 
c.Label.String = 'Relative Assembly Power'; 
c.Label.FontSize = 11; 
title('BIBLIS Assembly-Averaged Power Distribution'); 
xlabel('Assembly Column'); 
ylabel('Assembly Row'); 
xticks(1:17); yticks(1:17); 

% Overlay the exact numerical values onto the heatmap 
for row = 1:17 
    for col = 1:17 
        if active_assemblies(row, col) 
            % Determine text color based on background intensity for contrast 
            if relative_power(row, col) > 1.2 || relative_power(row, col) < 0.3 
                text_color = 'w'; 
            else 
                text_color = 'k'; 
            end 
             
            text(col, row, sprintf('%.3f', relative_power(row, col)), ... 
                'HorizontalAlignment', 'center', ... 
                'VerticalAlignment', 'middle', ... 
                'FontSize', 7, 'Color', text_color, 'FontWeight', 'bold'); 
        end 
    end 
end 

%% ======================================================================== 
%  9. COMPARISON PLOT ANALOGOUS TO ARTICLE FIGURE 6 
% ========================================================================= 
fprintf('\nGenerating Figure 6 analogue...\n'); 

N_compare   = N_refine(1); % Compares Coarsest (N=1) to Finest (N=6)
test_power = power_history_17x17(:, :, 1); 
ref_power  = power_history_17x17(:, :, 2); 

figure('Name', sprintf('Figure 6 Analogue: N=%d vs Reference', N_compare), ... 
       'Color', 'w', 'Position', [150, 150, 900, 850]); 
hold on; 
axis equal off; 
title(sprintf('BIBLIS Assembly-Averaged Power Distribution\nTop: FDM N=%d  |  Bottom: FDM Reference (N=%d)', ... 
      N_compare, N_refine(end)), 'FontSize', 12, 'FontWeight', 'bold'); 

for row = 1:17 
    for col = 1:17 
        if active_assemblies(row, col) 
            x_c = col; 
            y_c = 18 - row;  
             
            x0 = x_c - 0.5; x1 = x_c + 0.5; 
            y0 = y_c - 0.5; y1 = y_c + 0.5; 
            patch('XData', [x0 x1 x1 x0], 'YData', [y0 y0 y1 y1], ... 
                  'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8], 'LineWidth', 0.5); 
             
            text(x_c, y_c + 0.20, sprintf('%.4f', test_power(row, col)), ... 
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ... 
                'FontSize', 7, 'Color', 'b', 'FontWeight', 'bold'); 
             
            text(x_c, y_c - 0.20, sprintf('%.4f', ref_power(row, col)), ... 
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ... 
                'FontSize', 7, 'Color', 'k'); 
        end 
    end 
end 
xlim([0.5, 17.5]); ylim([0.5, 17.5]); 

%% ======================================================================== 
%  10. NORMALIZED LOWER-RIGHT QUADRANT COMPARISON (ARTICLE VS. MATLAB) 
% ========================================================================= 
fprintf('\nGenerating Normalized Lower-Right Quadrant Comparison...\n'); 
anm_full = NaN(17, 17); 
anm_data = { 
    [1.0909, 1.1015, 1.2426, 1.2205, 1.0886, 0.9820, 1.0947, 1.0145], ... 
    [1.1015, 1.1174, 1.1341, 1.2235, 1.0677, 1.0318, 1.0717, 0.9705], ... 
    [1.2426, 1.1341, 1.1221, 1.1050, 1.1199, 0.9236, 0.9309, 0.8245], ... 
    [1.2205, 1.2235, 1.1050, 1.1607, 1.0387, 0.9501, 0.7651, 0.5456], ... 
    [1.0886, 1.0677, 1.1199, 1.0387, 1.1226, 0.9931, 0.8746],         ... 
    [0.9820, 1.0318, 0.9236, 0.9501, 0.9931, 1.1997, 0.6845],         ... 
    [1.0947, 1.0718, 0.9309, 0.7651, 0.8746, 0.6845],                 ... 
    [1.0145, 0.9705, 0.8245, 0.5456]                                      
}; 
for i = 1:8 
    global_row = 8 + i;  
    row_vals = anm_data{i}; 
    anm_full(global_row, 9 : 8 + length(row_vals)) = row_vals; 
end 

anm_max = max(anm_full(:)); 
anm_norm = anm_full / anm_max; 

my_power = power_history_17x17(:, :, 2); % Use finest mesh
my_power(~active_assemblies) = NaN; 
my_norm = my_power / max(my_power(:)); 

figure('Name', 'Lower-Right Quadrant Comparison', 'Color', 'w', ... 
       'Position', [100, 150, 1200, 550]); 
cmap = jet(256); 

subplot(1, 2, 1);  
hold on; axis equal off; 
title('Article ANM 4x4 (Normalized)', 'FontSize', 14, 'FontWeight', 'bold'); 
for global_row = 9:17 
    for global_col = 9:17 
        if active_assemblies(global_row, global_col) && ~isnan(anm_norm(global_row, global_col)) 
            val = anm_norm(global_row, global_col); 
            color_idx = max(1, min(256, round(val * 255) + 1)); 
             
            x_center = global_col - 8; 
            y_center = 10 - (global_row - 8);  
             
            x0 = x_center - 0.5; x1 = x_center + 0.5; 
            y0 = y_center - 0.5; y1 = y_center + 0.5; 
             
            patch('XData', [x0 x1 x1 x0], 'YData', [y0 y0 y1 y1], ... 
                  'FaceColor', cmap(color_idx, :), 'EdgeColor', 'k'); 
                   
            text(x_center, y_center, sprintf('%.3f', val), ... 
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ... 
                'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold', ... 
                'BackgroundColor', [1 1 1 0.7]); 
        end 
    end 
end 
xlim([0.5, 9.5]); ylim([0.5, 9.5]); 

subplot(1, 2, 2);  
hold on; axis equal off; 
title(sprintf('Your MATLAB FDM N=%d (Normalized)', N_refine(end)), 'FontSize', 14, 'FontWeight', 'bold'); 
for global_row = 9:17 
    for global_col = 9:17 
        if active_assemblies(global_row, global_col) && ~isnan(my_norm(global_row, global_col)) 
            val = my_norm(global_row, global_col); 
            color_idx = max(1, min(256, round(val * 255) + 1)); 
             
            x_center = global_col - 8; 
            y_center = 10 - (global_row - 8);  
             
            x0 = x_center - 0.5; x1 = x_center + 0.5; 
            y0 = y_center - 0.5; y1 = y_center + 0.5; 
             
            patch('XData', [x0 x1 x1 x0], 'YData', [y0 y0 y1 y1], ... 
                  'FaceColor', cmap(color_idx, :), 'EdgeColor', 'k'); 
                   
            text(x_center, y_center, sprintf('%.3f', val), ... 
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ... 
                'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold', ... 
                'BackgroundColor', [1 1 1 0.7]); 
        end 
    end 
end 
xlim([0.5, 9.5]); ylim([0.5, 9.5]); 

colormap jet; caxis([0, 1]); 
cb = colorbar('Position', [0.93 0.15 0.015 0.7]); 
cb.Label.String = 'Normalized Power Distribution'; 
cb.FontSize = 11; 

%% ======================================================================== 
%  11. RELATIVE ERROR MAP (ANALOGOUS TO FEMFFUSION FIGURE 5) 
% ========================================================================= 
fprintf('\nGenerating Figure 5 Analogue (Relative Error %%)...\n'); 
P_test = power_history_17x17(:, :, 1);  
P_ref  = power_history_17x17(:, :, 2); 
error_pct = 100 * abs(P_test - P_ref) ./ (P_ref + eps); 
error_pct(~active_assemblies) = NaN; 

figure('Name', 'Relative Power Error (%)', 'Color', 'w', ... 
       'Position', [200, 200, 700, 600]); 
imagesc(error_pct, 'AlphaData', ~isnan(error_pct)); 
axis image;           
colormap jet; 
max_err = max(error_pct(active_assemblies)); 
clim([0, max_err]); 
c = colorbar; 
c.Label.String = 'Relative Error (%)'; 
c.Label.FontSize = 11; 
title(sprintf('Relative Error in Assembly Power\nFDM N=%d vs Reference N=%d', ... 
      N_refine(1), N_refine(end)), 'FontWeight', 'bold'); 
xlabel('Assembly Column'); 
ylabel('Assembly Row'); 
xticks(1:17); yticks(1:17);