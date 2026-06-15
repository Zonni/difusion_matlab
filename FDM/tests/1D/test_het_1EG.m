function test_het_1EG()
    load("het_1EG_exact_results.mat","size_vec","k_analytic", "phi_analytic_array","k_eff_num_vec");

    % PROBLEM PARAMETERS (1G Heterogeneous)
    L = 350;    
    x_int1 = 25;   % Coordinate of Interface 1 (Reflector -> Core)
    x_int2 = 325;  % Coordinate of Interface 2 (Core -> Reflector)
    
    % Physical properties
    % Index 1 = Reflector, Index 2 = Homogeneous Fuel
    D_lib          = [1.446; 0.776];
    sigma_a_lib    = [0.0077; 0.0244];
    nu_sigma_f_lib = [0; 0.0260];
    chi_lib        = [1; 1];
    sigma_s_lib = zeros(2, 1, 1);

    for i = 1:length(size_vec)  
        N = size_vec(i);
        
        region_lengths = [x_int1; x_int2 - x_int1; L - x_int2];
        region_materials = [1; 2; 1];
        cells_per_region = ceil(N / 3) * ones(3,1);
        
        mesh = Mesh_1D_FDM(region_lengths, cells_per_region);
        
        materials = NuclearData_1D(region_materials, D_lib, sigma_a_lib, nu_sigma_f_lib, chi_lib, sigma_s_lib);
        solver = Solver_1D_FDM(mesh, materials);
        solver = solver.assembleMatrices().solveEigenvalues(1);
        
        % Extract the computed k_eff
        k_eff_num = solver.keff(1);

        % Extract the computed flux
        flux_num = solver.phi(:, 1, 1);

        % Solver validation
        assert(abs(k_analytic - k_eff_num) < 1e-4); % Pretty bad comparison not to trigger errors
        
        assert(rmse(phi_analytic_array{i}, flux_num) < 1e-3);
    end
end