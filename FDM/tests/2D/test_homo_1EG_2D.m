function test_homo_1EG_2D()
%TEST_HOMO_1EG_2D  Convergence test for a 1-group homogeneous 2D reactor.
%
%   Loads the analytic k-eff from the pre-generated .mat file produced by
%   homo_1EG_2D_generate.m and asserts that the numerical solver converges
%   within the tolerance defined below for each mesh size.
%
%   Analytic reference:
%       k = nu_sigf / (sigma_a + 2*D*(pi/L)^2)
%

    addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', '2D'));

    %% Load data
    mat_path = fullfile(fileparts(mfilename('fullpath')), 'data', ...
                        'homo_1EG_2D_exact_results.mat');

    if ~isfile(mat_path)
        error(['Reference file not found: %s\n' ...
               'Run homo_1EG_2D_generate.m first.'], mat_path);
    end

    ref = load(mat_path, 'k_analytic');
    k_analytic = ref.k_analytic;

    %% Problem parameters
    L = 350;
    D_val = 0.776;
    SigA_val = 0.0244;
    NuSigF_val = 0.0260;

    D_lib = D_val;
    sigma_a_lib = SigA_val;
    nu_sigf_lib = NuSigF_val;
    chi_lib = 1.0;
    sigma_s_lib = 0.0;

    %% Convergence loop
    % Tolerance in pcm per mesh size.
    N_array = [40, 80, 120, 160];
    tol_array = [10, 10, 10, 10]; % pcm

    fprintf('Analytic k-eff: %.10f\n\n', k_analytic);
    fprintf('%-8s  %-12s  %-12s  %-12s  %-10s\n', ...
            'N', 'k_analytic', 'k_numeric', 'err (pcm)', 'tol (pcm)');
    fprintf('%s\n', repmat('-', 1, 60));

    for i = 1:length(N_array)
        N = N_array(i);
        tol = tol_array(i);

        region_lengths_x = (L / N) * ones(N, 1);
        cells_per_region_x = ones(N, 1);
    
        region_lengths_y = (L / N) * ones(N, 1);
        cells_per_region_y = ones(N, 1);

        mesh = Mesh_2D_FDM(region_lengths_x, cells_per_region_x, ...
                                   region_lengths_y, cells_per_region_y);

        region_materials = ones(N, N); % (Ny x Nx) — all active, 1 material

        nuclear_data = NuclearData_2D(region_materials, D_lib, sigma_a_lib, ...
                                      nu_sigf_lib, chi_lib, sigma_s_lib);

        solver = Solver_2D_FDM(mesh, nuclear_data);
        solver = solver.assembleMatrices();
        solver = solver.solveEigenvalues(1);

        k_num   = solver.keff(1);
        err_pcm = abs(k_analytic - k_num) * 1e5;

        fprintf('%-8d  %-12.5f  %-12.5f  %-12.1f  %-10d\n', ...
                N, k_analytic, k_num, err_pcm, tol);

        assert(err_pcm < tol, ...
            sprintf('N=%d: error %.1f pcm exceeds tolerance %d pcm.', ...
                    N, err_pcm, tol));
    end

    fprintf('\nAll assertions passed.\n');
end