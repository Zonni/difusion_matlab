function test_homo_2EG()
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')); 
    % --- Material properties (Homogeneous, 2 Energy Groups) ---
    % Index 1 = Fast Group (Group 1), Index 2 = Thermal Group (Group 2)
    D          = [1.446; 0.776];
    Sigma_a    = [0.0077; 0.0244];
    nu_Sigma_f = [0.0; 0.0260];
    
    % Group 1 to Group 2 down-scattering 
    % MADE UP VALUE
    Sigma_s12  = 0.015; 
    
    % Removal Cross Sections
    Sigma_r1 = Sigma_a(1) + Sigma_s12;
    Sigma_r2 = Sigma_a(2);
    
    % Cell information
    Delta_x = [25; 150; 150; 25];
    
    % Coupling Coefficients (d_g,i)
    % Since it is homogeneous, D_left = D_right = D(g).
    % Formula simplifies to: 2*D*D / (Dx_i*D + Dx_{i+1}*D) = 2*D / (Dx_i + Dx_{i+1})
    d1 = zeros(3,1);
    d2 = zeros(3,1);
    for i = 1:3
        d1(i) = 2 * D(1) / (Delta_x(i) + Delta_x(i+1));
        d2(i) = 2 * D(2) / (Delta_x(i) + Delta_x(i+1));
    end
    
    % --- Vacuum Boundary conditions (L_g,i) ---
    L1_left  = 2 * D(1) / Delta_x(1);
    L1_right = 2 * D(1) / Delta_x(4);
    
    L2_left  = 2 * D(2) / Delta_x(1);
    L2_right = 2 * D(2) / Delta_x(4);
    
    % --- Matrix Assembly (4x4 Blocks) ---
    A11 = zeros(4,4); % Destruction/Leakage Group 1
    A22 = zeros(4,4); % Destruction/Leakage Group 2
    A21 = zeros(4,4); % Scattering G1 -> G2
    B11 = zeros(4,4); % Fission from G1
    B12 = zeros(4,4); % Fission from G2
    
    % Block A11
    A11(1,1) = d1(1) + L1_left + Delta_x(1)*Sigma_r1;
    A11(1,2) = -d1(1);
    
    A11(2,1) = -d1(1);
    A11(2,2) = d1(1) + d1(2) + Delta_x(2)*Sigma_r1;
    A11(2,3) = -d1(2);
    
    A11(3,2) = -d1(2);
    A11(3,3) = d1(2) + d1(3) + Delta_x(3)*Sigma_r1;
    A11(3,4) = -d1(3);
    
    A11(4,3) = -d1(3);
    A11(4,4) = d1(3) + L1_right + Delta_x(4)*Sigma_r1;
    
    % Block A22
    A22(1,1) = d2(1) + L2_left + Delta_x(1)*Sigma_r2;
    A22(1,2) = -d2(1);
    
    A22(2,1) = -d2(1);
    A22(2,2) = d2(1) + d2(2) + Delta_x(2)*Sigma_r2;
    A22(2,3) = -d2(2);
    
    A22(3,2) = -d2(2);
    A22(3,3) = d2(2) + d2(3) + Delta_x(3)*Sigma_r2;
    A22(3,4) = -d2(3);
    
    A22(4,3) = -d2(3);
    A22(4,4) = d2(3) + L2_right + Delta_x(4)*Sigma_r2;
    
    % Blocks A21, B11, B12 (Diagonal matrices)
    for i = 1:4
        A21(i,i) = Delta_x(i) * Sigma_s12;
        B11(i,i) = Delta_x(i) * nu_Sigma_f(1);
        B12(i,i) = Delta_x(i) * nu_Sigma_f(2);
    end
    
    % --- Assemble Full 8x8 Matrices ---
    A = [ A11, zeros(4,4);
         -A21, A22 ];
         
    F = [ B11, B12;
          zeros(4,4), zeros(4,4) ];
          
    % --- Solve the Eigenvalue Problem ---
    % The system is A * Phi = (1/k) * F * Phi
    [V, D_eig] = eig(F, A);
    
    % Extract the eigenvalues
    eigenvalues = diag(D_eig);
    
    % Handle mathematical artifacts (eig might return NaNs for zero-blocks)
    valid_idx = isfinite(eigenvalues);
    eigenvalues(~valid_idx) = 0;
    
    % Find the dominant eigenvalue (k_eff) and corresponding eigenvector (flux)
    [k_eff_exact, idx] = max(abs(eigenvalues));
    phi_exact = V(:, idx);
    
    % Normalize the flux so the maximum value is 1
    phi_exact = phi_exact / max(abs(phi_exact));
    
    % Ensure positivity
    if sum(phi_exact) < 0
        phi_exact = -phi_exact;
    end
    
    % --- Display Results ---
    fprintf('Effective Multiplication Factor (k_eff) = %.5f\n\n', k_eff_exact);
    disp('Normalized Neutron Flux array:');
    disp('Node      Group 1 (Fast)      Group 2 (Thermal)');
    disp('--------------------------------------------------');
    for i = 1:4
        fprintf('%d         %10.5e         %10.5e\n', i, phi_exact(i), phi_exact(i+4));
    end
    
    
    %% 1D 2EG Homogeneous Problem (4 nodes)
    
    % MESH DATA
    % 4 regions, each representing a single node to match the 4-node problem
    region_lengths = [25; 150; 150; 25];
    cells_per_region = [1; 1; 1; 1];
    
    % NUCLEAR DATA
    % Homogeneous reactor -> all regions use Material 1
    region_materials = [1; 1; 1; 1]; 
    
    % Library dimensions convention: (num_materials x num_groups)
    % Column 1 = Fast Group (Group 1), Column 2 = Thermal Group (Group 2)
    D_lib       = [1.446, 0.776];  % in cm
    sigma_a_lib = [0.0077, 0.0244]; % in 1/cm
    nu_sigf_lib = [0.0000, 0.0260]; % Fission only in thermal group; 1/cm
    
    % Fission spectrum (Chi): all fission neutrons are born in the fast group
    chi_lib     = [1.0, 0.0];      
    
    % Scattering matrix: sigma_s_lib(material, g_to, g_from)
    % We have 1 material, 2 groups -> 1x2x2 matrix
    sigma_s_lib = zeros(1, 2, 2);
    
    % Down-scattering from Group 1 (fast) to Group 2 (thermal)
    sigma_s_lib(1, 2, 1) = 0.015;  % in 1/cm
    
    %% Material grid and mesh initialization
    mesh = Mesh_1D_FDM(region_lengths, cells_per_region);
    mesh.displayMesh();
    
    materials = NuclearData_1D(region_materials, D_lib, sigma_a_lib, nu_sigf_lib, chi_lib, sigma_s_lib);
    materials.displayMaterials();
    
    %% Solver initialization
    solver = Solver_1D_FDM(mesh, materials);
    
    % Assembling the system and solving the eigenvalue problem
    % Requesting 1 eigenvalue (fundamental mode)
    solver = solver.assembleMatrices().solveEigenvalues(1);
    
    k_eff_num = solver.keff(1);
    phi_num = solver.phi;
    A_num = solver.A;
    phi_num_internal = solver.phi(2:end-1, :, 1); 

    phi_exact_matrix = reshape(phi_exact, [4, 2]);
    
    max_error = max(abs(phi_exact_matrix - phi_num_internal), [], 'all');
    
    assert(max_error < 1e-10);
    assert(abs(k_eff_exact - k_eff_num) < 1e-12)
    assert(all(all(abs(A - A_num) < 1e-12)))    
end