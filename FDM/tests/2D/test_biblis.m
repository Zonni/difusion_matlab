function test_biblis()
%TEST_BIBLIS  Test for the Biblis benchmark.
%
%   Loads the analytic k-eff from the pre-generated .mat file produced by
%   biblis.m and asserts that the solver converges within the tolerance
%   defined below for each mesh size.
%
%   Benchmark reference:
%       k = 1.025110
%
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', '2D'));
    
    %% Load data
    mat_path = fullfile(fileparts(mfilename('fullpath')), 'data', 'biblis.mat');
    if ~isfile(mat_path)
        error(['Reference file not found: %s\n' ...
               'Run biblis.m first.'], mat_path);
    end
    ref = load(mat_path, 'k_ref');
    k_ref = ref.k_ref;
    
    %% ======================================================================== 
    %  NUCLEAR DATA  (fixed across all refinement levels) 
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
    %  ASSEMBLY-LEVEL MATERIAL GRID  (17 x 17) 
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
    %  ASSERTION LOOP 
    % ========================================================================= 
    N_refine = [1, 2, 3, 4, 5, 6]; % sub-cells per assembly side 
    n_levels = length(N_refine); 
    
    keff_vec = zeros(n_levels, 1); 
    dof_vec  = zeros(n_levels, 1); 
    err_pcm  = zeros(n_levels, 1);
    
    % Define specific error tolerances (pcm) for each refinement level
    tol_vec  = [1600, 350, 150, 80, 50, 30]; 
    
    fprintf('%-6s  %-12s  %-10s  %-12s  %-12s  %-10s\n', ...
            'N', 'Grid', 'DOF', 'k_numeric', 'err (pcm)', 'tol (pcm)');
    fprintf('%s\n', repmat('-', 1, 70));
    
    for refinment = 1:n_levels 
        N = N_refine(refinment); 
         
        % Expand assembly grid 
        region_materials = kron(assembly_mat, ones(N, N)); 
         
        % Mesh 
        region_lengths = assembly_size * ones(17, 1); 
        cells_per_region = N * ones(17, 1); 
        mesh = Mesh_2D_FDM(region_lengths, cells_per_region, ... 
                           region_lengths, cells_per_region); 
                        
        nuclear_data = NuclearData_2D(region_materials, D_lib, sigma_a_lib, ... 
                                      nu_sigf_lib, chi_lib, sigma_s_lib); 
                                   
        solver = Solver_2D_FDM(mesh, nuclear_data, 'vacuum'); 
        solver = solver.assembleMatrices(); 
         
        solver = solver.solveEigenvalues(1); 
         
        % Store current metrics
        k_num = solver.keff(1);
        keff_vec(refinment) = k_num; 
        dof_vec(refinment) = 2 * numel(solver.active_cells);    
        
        current_err = abs(k_num - k_ref) * 1e5; 
        err_pcm(refinment) = current_err;
        
        % Check against the tightening tolerance
        tol = tol_vec(refinment);
         
        fprintf('%-6d  %-12s  %-10d  %-12.6f  %-12.1f  %-10d\n', ... 
                N, sprintf('%dx%d', 17*N, 17*N), dof_vec(refinment), ...
                k_num, current_err, tol);
        
        assert(current_err < tol, ...
            sprintf('N=%d: error %.1f pcm exceeds tolerance %d pcm.', ...
                    N, current_err, tol));
    end 
    
    fprintf('\nAll assertions passed successfully!\n');
end