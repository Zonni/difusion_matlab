classdef NuclearData_2D_FDM
    %NUCLEARDATA_2D Nuclear material properties for 2D multigroup diffusion.

    properties
        num_materials    (1,1) double
        num_groups       (1,1) double
        region_materials (:,:) double   % (Ny x Nx)  0 = outside reactor
        D_lib            (:,:) double   % (material, group)
        sigma_a_lib      (:,:) double   % (material, group)
        nu_sigf_lib      (:,:) double   % (material, group)
        chi_lib          (:,:) double   % (material, group)
        sigma_s_lib      (:,:,:) double % (material, g_from, g_to)
    end

    methods
        function obj = NuclearData_2D_FDM(region_materials, D_lib, sigma_a_lib, ...
                                      nu_sigf_lib, chi_lib, sigma_s_lib)
            [n_mat, n_grp] = size(D_lib);

            if ~isequal(size(sigma_a_lib), [n_mat, n_grp]) || ...
               ~isequal(size(nu_sigf_lib), [n_mat, n_grp]) || ...
               ~isequal(size(chi_lib),     [n_mat, n_grp])
                error("D_lib, sigma_a_lib, nu_sigf_lib, and chi_lib must all be (%d x %d).", n_mat, n_grp);
            end

            if size(sigma_s_lib, 1) ~= n_mat || ...
               size(sigma_s_lib, 2) ~= n_grp || ...
               size(sigma_s_lib, 3) ~= n_grp
                error("sigma_s_lib must be (%d x %d x %d).", n_mat, n_grp, n_grp);
            end

            max_idx = max(region_materials(:));
            if max_idx > n_mat
                error("region_materials references material %d but only %d defined.", max_idx, n_mat);
            end
            if any(region_materials(:) < 0)
                error("region_materials must be non-negative (0 = outside reactor).");
            end

            obj.num_materials    = n_mat;
            obj.num_groups       = n_grp;
            obj.region_materials = region_materials;
            obj.D_lib            = D_lib;
            obj.sigma_a_lib      = sigma_a_lib;
            obj.nu_sigf_lib      = nu_sigf_lib;
            obj.chi_lib          = chi_lib;
            obj.sigma_s_lib      = reshape(sigma_s_lib, n_mat, n_grp, n_grp);
        end

        function sigma_r = getRemovalCrossSection(obj, mat_id, g)
            % sigma_r,g = sigma_a,g + sum_{h~=g} sigma_s,g->h
            g_others = setdiff(1:obj.num_groups, g);
            sigma_r  = obj.sigma_a_lib(mat_id, g) + ...
                       sum(obj.sigma_s_lib(mat_id, g, g_others), 3);
        end

        function displayMaterials(obj)
            fprintf('======================\n  MATERIALS DATA  \n======================\n\n');
            fprintf('Unique Materials : %d\n',   obj.num_materials);
            fprintf('Energy Groups    : %d\n\n', obj.num_groups);
            for m = 1:obj.num_materials
                fprintf('--- Material %d ---\n', m);
                fprintf('  D       : '); fprintf('%8.4f  ', obj.D_lib(m,:));       fprintf('\n');
                fprintf('  sigma_a : '); fprintf('%8.4f  ', obj.sigma_a_lib(m,:)); fprintf('\n');
                fprintf('  nu_sigf : '); fprintf('%8.4f  ', obj.nu_sigf_lib(m,:)); fprintf('\n');
                fprintf('  chi     : '); fprintf('%8.4f  ', obj.chi_lib(m,:));     fprintf('\n');
            end
            fprintf('\nMaterial grid (0 = outside reactor):\n');
            disp(obj.region_materials);
        end
    end
end