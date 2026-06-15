function test_homo_1EG()
    load("homo_1EG_exact_results.mat", "N_array", "k_analytic", "k_eff_num_vec");

    % PROBLEM PARAMETERS (1G Homogeneous)
    L = 350;
        
    % Physical properties (Homogeneous Fuel)
    D_val          = 0.776;   % cm
    sigma_a_val    = 0.0244;  % 1/cm
    nu_sigma_f_val = 0.0260;  % 1/cm
    sigma_s_val    = zeros(1, 1, 1);
    chi_val = 1;

    region_materials = 1;  % single region, single material
    
    for i = 1:length(N_array)
        N = N_array(i);
        mesh = Mesh_1D_FDM(L, N);
        materials = NuclearData_1D(region_materials, D_val, sigma_a_val, ...
                              nu_sigma_f_val, chi_val, sigma_s_val);
        solver = Solver_1D_FDM(mesh, materials);
        solver = solver.assembleMatrices().solveEigenvalues(1);

        k_eff_num = solver.keff(1);
        assert(abs(k_analytic - k_eff_num) < 1e-5);
    end
end
