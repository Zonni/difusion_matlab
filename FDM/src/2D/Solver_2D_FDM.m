classdef Solver_2D_FDM
    %SOLVER_2D_FDM Multigroup 2D neutron diffusion solver.
    %   Solves A*phi = (1/k)*B*phi using mesh-centred FDM (Hérbert).
    %
    %   Boundary conditions
    %   -------------------
    %   'zero_flux' : phi = 0 at the physical boundary face.
    %                 Boundary coefficient: 2*D*transverse / delta
    %   'vacuum'    : phi = 0 at extrapolated distance 2D beyond the face.
    %                 Boundary coefficient: 2*D*transverse / (delta + 4*D)

    properties
        mesh         Mesh_2D_FDM
        nuclear_data NuclearData_2D_FDM
        bc_type      char            % 'zero_flux' | 'vacuum'
        A            (:,:)   double
        B            (:,:)   double
        phi          (:,:,:) double  % (N_active, group, mode)
        keff         (:,1)   double
        power
        active_cells (:,1)   double
    end

    methods
        function obj = Solver_2D_FDM(mesh, nuclear_data, bc_type)
            %SOLVER_2D_FDM Construct the solver.
            %
            %   Inputs:
            %       mesh         - Mesh_2D_FDM object.
            %       nuclear_data - NuclearData_2D object.
            %       bc_type      - (Optional) 'vacuum' (default) or 'zero_flux'.

            [Ny_mat, Nx_mat] = size(nuclear_data.region_materials);
            if Nx_mat ~= mesh.num_cells_x || Ny_mat ~= mesh.num_cells_y
                error("region_materials is (%d x %d) but mesh has " + ...
                      "num_cells_y = %d, num_cells_x = %d.", ...
                      Ny_mat, Nx_mat, mesh.num_cells_y, mesh.num_cells_x);
            end

            obj.mesh         = mesh;
            obj.nuclear_data = nuclear_data;
            if nargin < 3 || isempty(bc_type)
                obj.bc_type = 'vacuum';
            else
                obj.bc_type = bc_type;
            end
        end

        % -----------------------------------------------------------------
        function obj = assembleMatrices(obj)
            %ASSEMBLEMATRICES Build the loss matrix A and fission matrix B.

            Nx      = obj.mesh.num_cells_x;
            Ny      = obj.mesh.num_cells_y;
            G       = obj.nuclear_data.num_groups;
            reg_mat = obj.nuclear_data.region_materials;   % (Ny x Nx)
            vacuum  = strcmp(obj.bc_type, 'vacuum');

            % Full-to-active index map
            obj.active_cells = find(reg_mat(:) ~= 0);
            N_active         = numel(obj.active_cells);
            active_idx       = zeros(Nx * Ny, 1);
            active_idx(obj.active_cells) = 1:N_active;

            A = sparse(G * N_active, G * N_active);
            B = sparse(G * N_active, G * N_active);

            for g = 1:G
                g_off = (g - 1) * N_active;

                for k_active = 1:N_active
                    k_full          = obj.active_cells(k_active);
                    m_C             = reg_mat(k_full);
                    D_C             = obj.nuclear_data.D_lib(m_C, g);
                    [row_C, col_C]  = ind2sub([Ny, Nx], k_full);
                    dx_C            = obj.mesh.delta_x(col_C);
                    dy_C            = obj.mesh.delta_y(row_C);

                    diag_leak = 0;

                    % --- Left face (x-) ---
                    left_is_bnd = (col_C == 1) || ...
                        (reg_mat(sub2ind([Ny,Nx], row_C, col_C-1)) == 0);
                    if left_is_bnd
                        diag_leak = diag_leak + obj.bnd_coeff(D_C, dy_C, dx_C, vacuum);
                    else
                        k_L  = active_idx(sub2ind([Ny,Nx], row_C, col_C-1));
                        D_L  = obj.nuclear_data.D_lib(reg_mat(sub2ind([Ny,Nx], row_C, col_C-1)), g);
                        dx_L = obj.mesh.delta_x(col_C - 1);
                        aL   = obj.int_coeff(D_C, D_L, dy_C, dx_C, dx_L);
                        A(g_off + k_active, g_off + k_L) = -aL;
                        diag_leak = diag_leak + aL;
                    end

                    % --- Right face (x+) ---
                    right_is_bnd = (col_C == Nx) || ...
                        (reg_mat(sub2ind([Ny,Nx], row_C, col_C+1)) == 0);
                    if right_is_bnd
                        diag_leak = diag_leak + obj.bnd_coeff(D_C, dy_C, dx_C, vacuum);
                    else
                        k_R  = active_idx(sub2ind([Ny,Nx], row_C, col_C+1));
                        D_R  = obj.nuclear_data.D_lib(reg_mat(sub2ind([Ny,Nx], row_C, col_C+1)), g);
                        dx_R = obj.mesh.delta_x(col_C + 1);
                        aR   = obj.int_coeff(D_C, D_R, dy_C, dx_C, dx_R);
                        A(g_off + k_active, g_off + k_R) = -aR;
                        diag_leak = diag_leak + aR;
                    end

                    % --- Top face (y-) ---
                    top_is_bnd = (row_C == 1) || ...
                        (reg_mat(sub2ind([Ny,Nx], row_C-1, col_C)) == 0);
                    if top_is_bnd
                        diag_leak = diag_leak + obj.bnd_coeff(D_C, dx_C, dy_C, vacuum);
                    else
                        k_T  = active_idx(sub2ind([Ny,Nx], row_C-1, col_C));
                        D_T  = obj.nuclear_data.D_lib(reg_mat(sub2ind([Ny,Nx], row_C-1, col_C)), g);
                        dy_T = obj.mesh.delta_y(row_C - 1);
                        aT   = obj.int_coeff(D_C, D_T, dx_C, dy_C, dy_T);
                        A(g_off + k_active, g_off + k_T) = -aT;
                        diag_leak = diag_leak + aT;
                    end

                    % --- Bottom face (y+) ---
                    bot_is_bnd = (row_C == Ny) || ...
                        (reg_mat(sub2ind([Ny,Nx], row_C+1, col_C)) == 0);
                    if bot_is_bnd
                        diag_leak = diag_leak + obj.bnd_coeff(D_C, dx_C, dy_C, vacuum);
                    else
                        k_B  = active_idx(sub2ind([Ny,Nx], row_C+1, col_C));
                        D_B  = obj.nuclear_data.D_lib(reg_mat(sub2ind([Ny,Nx], row_C+1, col_C)), g);
                        dy_B = obj.mesh.delta_y(row_C + 1);
                        aB   = obj.int_coeff(D_C, D_B, dx_C, dy_C, dy_B);
                        A(g_off + k_active, g_off + k_B) = -aB;
                        diag_leak = diag_leak + aB;
                    end

                    % --- Diagonal: leakage + removal ---
                    sigma_r = obj.nuclear_data.getRemovalCrossSection(m_C, g);
                    A(g_off + k_active, g_off + k_active) = diag_leak + sigma_r * dx_C * dy_C;

                    % --- In-scatter (off-diagonal A blocks) ---
                    for h = 1:G
                        if h ~= g
                            h_off    = (h - 1) * N_active;
                            sigma_hg = obj.nuclear_data.sigma_s_lib(m_C, h, g);
                            A(g_off + k_active, h_off + k_active) = ...
                                A(g_off + k_active, h_off + k_active) - sigma_hg * dx_C * dy_C;
                        end
                    end

                    % --- Fission matrix B ---
                    chi_g = obj.nuclear_data.chi_lib(m_C, g);
                    for h = 1:G
                        h_off = (h - 1) * N_active;
                        B(g_off + k_active, h_off + k_active) = ...
                            chi_g * obj.nuclear_data.nu_sigf_lib(m_C, h) * dx_C * dy_C;
                    end
                end
            end

            obj.A = A;
            obj.B = B;
        end

        % -----------------------------------------------------------------
        function obj = solveEigenvalues(obj, nm)
            %SOLVEEIGENVALUES Solve B*phi = k*A*phi for nm modes.
            [V, D_eig] = eigs(obj.B, obj.A, nm, 'largestabs');
            obj.keff   = diag(D_eig);

            V_norm = V ./ max(abs(V));
            for m = 1:nm
                if sum(V_norm(:, m)) < 0
                    V_norm(:, m) = -V_norm(:, m);
                end
            end

            N_active = numel(obj.active_cells);
            G        = obj.nuclear_data.num_groups;
            obj.phi  = zeros(N_active, G, nm);
            for g = 1:G
                idx_g = (g-1)*N_active + 1 : g*N_active;
                obj.phi(:, g, :) = V_norm(idx_g, :);
            end
        end

        % -----------------------------------------------------------------
        function displayProblem(obj)
            %DISPLAYPROBLEM Print the effective multiplication factors.
            fprintf('=================\n  PROBLEM DATA  \n=================\n\n');
            fprintf('Boundary condition: %s\n\n', obj.bc_type);
            fprintf('K-eff (Fundamental): %.5f\n', obj.keff(1));
            if length(obj.keff) > 1
                fprintf('K-eff (Harmonics):   ');
                fprintf('%.5f ', obj.keff(2:end));
                fprintf('\n');
            end
        end

        % -----------------------------------------------------------------
        function plotPhi(obj, mode_idx)
            %PLOTPHI Plot the 2D flux for every energy group.
            if nargin < 2
                mode_idx = 1;
            end
            Nx = obj.mesh.num_cells_x;
            Ny = obj.mesh.num_cells_y;
            G  = obj.nuclear_data.num_groups;
            x  = obj.mesh.cell_centers_x;
            y  = obj.mesh.cell_centers_y;
            
            for g = 1:G
                flux_full                   = NaN(Ny, Nx);
                flux_full(obj.active_cells) = obj.phi(:, g, mode_idx);
                figure('Name', sprintf('Flux — Group %d, Mode %d', g, mode_idx), ...
                       'Color', 'w');
                colormap jet
                imagesc(x, y, flux_full);
                set(gca, 'YDir', 'normal');
                colorbar;
                xlabel('X (cm)'); ylabel('Y (cm)');
                title(sprintf('Flux — Group %d | Mode %d | k_{eff} = %.5f', ...
                              g, mode_idx, obj.keff(mode_idx)));
                axis equal tight;
            end
        end
    
        function obj = computePower(obj, nm)
            % COMPUTEPOWER Computes the raw (unnormalized) spatial power distribution.
            if nargin < 2 || isempty(nm)
                nm = 1; 
            end
            
            Nx = obj.mesh.num_cells_x;
            Ny = obj.mesh.num_cells_y;
            G  = obj.nuclear_data.num_groups;
            num_modes = length(nm);
            
            obj.power = NaN(Ny, Nx, num_modes);
            m_vec = obj.nuclear_data.region_materials(obj.active_cells);
            cross_section_matrix = obj.nuclear_data.nu_sigf_lib(m_vec, 1:G); 
            
            for i = 1:num_modes
                n = nm(i); 
                phi_mode = abs(obj.phi(:, :, n));
                
                % 1. Calculate raw local power (P = sum(Sigma_f * Phi))
                power_active = sum(cross_section_matrix .* phi_mode, 2);
                
                % 2. Map the RAW power back to the 2D grid
                power_2D = NaN(Ny, Nx);
                power_2D(obj.active_cells) = power_active;
                
                obj.power(:, :, i) = power_2D;
            end
        end

        function plotPower(obj, nm)
            %PLOTPOWER Plot the normalized 2D power distribution for a specific mode
            if nargin < 2
                nm = 1;
            end
            
            pwr_map = obj.power(:,:,nm);
            
            figure('Name', sprintf('Power - Mode %d', nm), 'Color', 'w');
            imagesc(obj.mesh.cell_centers_x, obj.mesh.cell_centers_y, pwr_map, ...
                   'AlphaData', ~isnan(pwr_map));
            set(gca, 'YDir', 'normal', 'Color', [0.9 0.9 0.9]); 
            colormap jet;
            colorbar; 
            axis equal tight;
            xlabel('X (cm)'); ylabel('Y (cm)');
            title(sprintf('Neutron Power Distribution — Mode %d', nm));
        end    
    end
    % -----------------------------------------------------------------
    methods (Access = private, Static)
        function a = bnd_coeff(D, transverse, delta, vacuum)
            %BND_COEFF Boundary face leakage coefficient.
            %   vacuum = false : zero-flux at face    -> 2*D*T / delta
            %   vacuum = true  : zero-flux at 2D beyond face -> 2*D*T / (delta + 4*D)
            if vacuum
                a = 2 * D * transverse / (delta + 4 * D);
            else
                a = 2 * D * transverse / delta;
            end
        end

        function a = int_coeff(D_C, D_N, transverse, delta_C, delta_N)
            %INT_COEFF Interior face harmonic-mean coupling coefficient.
            a = (2 * D_C * D_N * transverse) / (delta_N * D_C + delta_C * D_N);
        end
    end
end