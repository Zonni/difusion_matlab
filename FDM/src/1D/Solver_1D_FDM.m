classdef Solver_1D_FDM
    %SOLVER_1D_FDM 1D multigroup neutron diffusion solver.
    %   Solves the steady-state eigenvalue problem using a mesh-centred 
    %   finite difference method approximation (Hérbert).
    properties
        mesh                  Mesh_1D_FDM
        nuclear_data          NuclearData_1D_FDM
        A                     (:,:)   double
        B                     (:,:)   double
        phi                   (:,:,:) double
        keff                  (:,1)   double
        full_mesh_coordinates (:,1) double
    end
    methods
        function obj = Solver_1D_FDM(mesh, nuclear_data)
            %SOLVER_1D_FDM Construct an instance of the 1D solver.
            %
            %   Inputs:
            %       mesh         - Mesh_1D_FDM object containing the spatial grid.
            %       nuclear_data - NuclearData_1D object containing the material properties.
            %
            %   Outputs:
            %       obj          - Initialized Solver_1D_FDM object.
            
            % Security size checks
            if length(mesh.region_lengths) ~= length(nuclear_data.region_materials)
                error("region_lengths and region_materials must have the same size.")
            end
            
            obj.mesh = mesh;
            obj.nuclear_data = nuclear_data;
            obj.full_mesh_coordinates = [mesh.region_boundaries(1); ...
                                         mesh.cell_centers; ...
                                         mesh.region_boundaries(end)];
        end

        function obj = assembleMatrices(obj)
            %ASSEMBLEMATRICES Builds the destruction (A) and production (B) matrices.
            %   Constructs the global block-sparse matrices representing leakage, 
            %   removal, scattering, and fission sources.
            %
            %   Inputs:
            %       obj - Solver_1D_FDM object.
            %
            %   Outputs:
            %       obj - Solver_1D_FDM object with populated A and B properties.
            
            % Assemble the matrices of the system

            num_cells = obj.mesh.num_cells;
            num_groups = obj.nuclear_data.num_groups;
            mat_ids = obj.nuclear_data.region_materials(obj.mesh.cell_region_idx);
            delta = obj.mesh.delta_x;
            total_dof = num_groups * num_cells;
            A = sparse(total_dof, total_dof);
            B = sparse(total_dof, total_dof);

            for g = 1:num_groups
                % Construct A_{g,g} block
                idx_g = (g-1)*num_cells + 1 : g*num_cells;

                D = obj.nuclear_data.D_lib(mat_ids, g);
                chi_g = obj.nuclear_data.chi_lib(mat_ids, g);
                sigma_r = obj.nuclear_data.getRemovalCrossSection(mat_ids, g);

                d_curr = D(1:end-1);
                d_next = D(2:end);
                dx_curr = delta(1:end-1);
                dx_next = delta(2:end);

                d = (2 .* d_curr .* d_next) ./ (dx_next .* d_curr + dx_curr .* d_next);

                l_left = (2 * D(1)) / delta(1);
                l_right = (2 * D(end)) / delta(end);

                diag_a = delta .* sigma_r + [d; 0] + [0; d];
                diag_a(1) = diag_a(1) + l_left;
                diag_a(end) = diag_a(end) + l_right;

                A(idx_g, idx_g) = spdiags(diag_a, 0, num_cells, num_cells) ...
                    + diag(sparse(-d),  1) ...
                    + diag(sparse(-d), -1);

                for h = 1:num_groups
                    % Construct B_{g,h} block
                    idx_h = (h-1)*num_cells + 1 : h*num_cells;
                    nu_sigf_h = obj.nuclear_data.nu_sigf_lib(mat_ids, h);
                    B(idx_g, idx_h) = spdiags(delta .* chi_g .* nu_sigf_h, 0, num_cells, num_cells);
                    
                    % Construct S_{g,h} block
                    if g ~= h
                        sigma_s_hg = reshape(obj.nuclear_data.sigma_s_lib(mat_ids, g, h), num_cells, 1);
                        A(idx_g, idx_h) = -spdiags(delta .* sigma_s_hg, 0, num_cells, num_cells);
                    end
                end
            end

            obj.A = A;
            obj.B = B;
        end
        function obj = solveEigenvalues(obj, nm)
            %SOLVEEIGENVALUES Solves the generalized algebraic eigenvalue problem.
            %
            %   Inputs:
            %       obj - Solver_1D_FDM object.
            %       nm  - Number of modes (eigenvalues) to compute.
            %
            %   Outputs:
            %       obj - Solver_1D_FDM object with populated keff and phi properties.

            [V, D_eig] = eigs(obj.B, obj.A, nm, 'largestabs');

            obj.keff = diag(D_eig);

            V_norm = V ./ max(abs(V));

            % Check if the sign is correct
            if sum(V_norm(:,1)) < 0
                V_norm(:,1) = -V_norm(:,1);
            end

            num_cells = obj.mesh.num_cells;
            num_groups = obj.nuclear_data.num_groups;
            obj.phi = zeros(num_cells + 2, num_groups, nm);

            for g = 1:num_groups
                idx_g = (g-1)*num_cells + 1 : g*num_cells;
                obj.phi(2:end-1, g, :) = V_norm(idx_g, :);
            end

        end
        function displayProblem(obj)
            %DISPLAYPROBLEM Prints the calculated effective multiplication factors.
            %
            %   Inputs:
            %       obj - Solver_1D_FDM object.
            %
            %   Outputs:
            %       None.
            fprintf('=================\n  PROBLEM DATA  \n=================\n\n');
            fprintf('K-eff (Fundamental): %.5f\n', obj.keff(1));
            if length(obj.keff) > 1
                fprintf('K-eff (Harmonics):   ');
                fprintf('%.5f ', obj.keff(2:end));
                fprintf('\n');
            end
        end

        function plotPhi(obj, mode_idx)
            %PLOTPHI Plots the spatial distribution of the neutron flux.
            %
            %   Inputs:
            %       obj      - Solver_1D_FDM object.
            %       mode_idx - (Optional) Index of the eigenmode to plot. Default is 1.
            %
            %   Outputs:
            %       None.
            
            if nargin < 2
                mode_idx = 1;
            end

            num_groups = obj.nuclear_data.num_groups;
            figure; hold on; grid on;
            colors = lines(num_groups);

            for g = 1:num_groups
                plot(obj.full_mesh_coordinates, obj.phi(:, g, mode_idx), ...
                    'LineWidth', 1.5, 'Color', colors(g,:), ...
                    'DisplayName', sprintf('Group %d', g));
            end
            
            hold off;
            xlabel('Position x (cm)');
            ylabel('Normalized Flux \Phi(x)');
            title(sprintf('Flux Profile — Mode %d  (k = %.5f)', mode_idx, obj.keff(mode_idx)));
            legend show;
        end
    end
end