classdef NuclearData_1D_FDM
    %NUCLEAR_DATA_1D Nuclear material properties for 1D multigroup diffusion equation.
    %   Stores and manages multigroup cross sections and diffusion coefficients.
    properties
        num_materials    (1,1) double
        num_groups       (1,1) double
        region_materials (:,1) double
        D_lib            (:,:) double       % (material, group)
        sigma_a_lib      (:,:) double       % (material, group)
        nu_sigf_lib      (:,:) double       % (material, group)
        chi_lib          (:,:) double       % (material, group)
        sigma_s_lib      (:,:,:) double     % (material, g_from, g_to)
    end

    methods
        function obj = NuclearData_1D_FDM(region_materials, D_lib, sigma_a_lib, nu_sigf_lib, chi_lib, sigma_s_lib)
            %NUCLEAR_DATA Construct an instance of the nuclear data manager.
            %
            %   Inputs:
            %       region_materials - (N x 1) Material index for each physical region.
            %       D_lib            - (M x G) Diffusion coefficients.
            %       sigma_a_lib      - (M x G) Macroscopic absorption cross sections.
            %       nu_sigf_lib      - (M x G) Nu times macroscopic fission cross sections.
            %       chi_lib          - (M x G) Fission spectrum.
            %       sigma_s_lib      - (M x G x G) Macroscopic scattering matrix (material, from_group, to_group).
            %
            %   Outputs:
            %       obj              - Initialized NuclearData_1D object.
                        
            % Security size checks
            if ~isequal(size(D_lib), size(sigma_a_lib), size(nu_sigf_lib), size(chi_lib))
                error("D_lib, sigma_a_lib, nu_sigf_lib, and chi_lib must all be (num_materials x num_groups).")
            end

            [n_mat, n_grp] = size(D_lib);

            if size(sigma_s_lib, 1) ~= n_mat || size(sigma_s_lib, 2) ~= n_grp || size(sigma_s_lib, 3) ~= n_grp
                error("sigma_s_lib must be (num_materials x num_groups x num_groups).")
            end
            
            obj.num_materials = n_mat;
            obj.num_groups = n_grp;
            obj.region_materials = region_materials;
            obj.D_lib = D_lib;
            obj.sigma_a_lib = sigma_a_lib;
            obj.nu_sigf_lib = nu_sigf_lib;
            obj.chi_lib = chi_lib;
            obj.sigma_s_lib = reshape(sigma_s_lib, n_mat, n_grp, n_grp);
        end

        function sigma_r = getRemovalCrossSection(obj, mat_ids, g)
            %GETREMOVALCROSSSECTION Calculates the macroscopic removal cross section.
            %
            %   Inputs:
            %       obj     - NuclearData_1D object.
            %       mat_ids - Array of material indices.
            %       g       - Current energy group index.
            %
            %   Outputs:
            %       sigma_r - Array of removal cross sections for the given materials and group.

            g_others = setdiff(1:obj.num_groups, g);
            
            sum_scatter_out = sum(obj.sigma_s_lib(mat_ids, g, g_others), 3);
            sigma_r = obj.sigma_a_lib(mat_ids, g) + sum_scatter_out;
            sigma_r = sigma_r(:);
        end
        
        function displayMaterials(obj)
            %DISPLAYMATERIALS Prints a summary of the macroscopic cross sections.
            %
            %   Inputs:
            %       obj - NuclearData_1D object.
            %
            %   Outputs:
            %       None.

            fprintf('======================\n  MATERIALS DATA  \n======================\n\n');
            fprintf('Unique Materials: %d\n',   obj.num_materials);
            fprintf('Energy Groups:    %d\n\n', obj.num_groups);
            for m = 1:obj.num_materials
                fprintf('--- Material %d ---\n', m);
                fprintf('  D:       '); fprintf('%.4f  ', obj.D_lib(m,:));       fprintf('\n');
                fprintf('  sigma_a: '); fprintf('%.4f  ', obj.sigma_a_lib(m,:)); fprintf('\n');
                fprintf('  nu_sigf: '); fprintf('%.4f  ', obj.nu_sigf_lib(m,:)); fprintf('\n');
                fprintf('  chi:     '); fprintf('%.4f  ', obj.chi_lib(m,:));     fprintf('\n');
            end
        end
    end
end