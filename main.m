clear; close all; clc;
data_juanjo;
% Save root_path to construc paths easily
root_path = fileparts(mfilename('fullpath'));
if method == "FEM"
        addpath(fullfile(root_path, "FEM"));

    % Load bank from data_bank (only physical data)
    bank_path = fullfile(root_path, 'data_bank');
    addpath(bank_path);
    run([char(bank), '.m']);
    rmpath(bank_path);

    % Collect bank variables into a struct
    bank_data = struct();
    vars_list = who();
    for i = 1:length(vars_list)
        var_name = vars_list{i};
        if ~ismember(var_name, {'bank_data', 'vars_list', 'root_path', 'bank_path', 'bank'})
            bank_data.(var_name) = eval(var_name);
        end
    end

    % Build params structure
    params = build_params(energy_groups, n_refinements, fe_degree, perturbation_function, ...
                          boundary_conditions, output_filename, output_flag, ...
                          transient, time_delta, time_end, ...
                          region_lengths, cells_per_region, region_materials);

    % Merge physical data from bank (with index transformation)
    params = merge_bank_data(params, bank_data);
    params.nombre = bank;

    % Build default time functions and apply perturbations
    params = time_functions(params);
    params = apply_perturbation(params, transient, perturbation_function, ...
                                material_changing, slope);

    % Message
    fprintf('\n========================================\n');
    fprintf('  SIMULATION\n');
    fprintf('========================================\n');
    fprintf('  Method      : %s\n', method);
    fprintf('  Solver      : %s\n', Solver);
    fprintf('  Bank       : %s\n', bank);
    fprintf('  Groups      : %d\n', energy_groups);
    fprintf('  Transient : %s\n', string(transient));
    fprintf('========================================\n\n');

    % Solver
    if transient
        [transitorio, t, y] = run_fem_solver(params, output_flag, Solver);
    else
        [transitorio, t, y] = run_stationary(params, output_flag);
    end

    % Save results
    save(params.output_filename, 't', 'y', 'transitorio');
    fprintf('Results saved in %s\n', params.output_filename);
    if output_flag
        fprintf('Graphs saved in results/figures/Resultados_*.eps\n');
    end
elseif method == "FDM"
    addpath(fullfile(root_path, "FDM", "src", "1D"));
    addpath(fullfile(root_path, "FDM", "src", "2D"));
    data_bank = fullfile(root_path, 'data_bank', bank);
    run(data_bank);
    
    if transient
        warning("FDM Branch: The FDM solver currently only supports steady-state eigenvalue " + ...
                "calculations. Proceeding with static eigenvalue execution.");
    end
    % Map boundary conditions to strings expected by FDM classes
    switch boundary_conditions
        case 0
            bc = 'zero_flux';
        case 1
            bc = 'vacuum';
        otherwise
            error("Invalid boundary_conditions = %d. FDM supports 0 (zero_flux) or 2 (vacuum).", ...
                  boundary_conditions);
    end
    % Dimensionality Cases
    switch dimension
        case 1
            fprintf("==> Initializing 1D Multi-Group FDM Discretization...\n");
                        
            % Instantiate geometric mesh object
            mesh = Mesh_1D_FDM(region_lengths, cells_per_region);
            
            % Instantiate physical data map
            nuclear_data = NuclearData_1D_FDM(region_materials, D_lib, sigma_a_lib, nu_sigf_lib, chi_lib, sigma_s_lib);
                                          
            % Instantiate 1D Solver environment
            solver = Solver_1D_FDM(mesh, nuclear_data);
        case 2
            fprintf("==> Initializing 2D Multi-Group FDM Discretization...\n");
            
            % 1. Automatically calculate the number of regions from the databank matrix
            n_regions_x = size(region_materials, 2); % 17
            n_regions_y = size(region_materials, 1); % 17
            
            % 2. Generate region lengths vectors assuming uniform assembly widths
            region_lengths_x = ones(n_regions_x, 1) * assembly_size;
            region_lengths_y = ones(n_regions_y, 1) * assembly_size;
            
            % 3. Calculate cells per region accounting for global grid refinements
            cells_x = ones(n_regions_x, 1) * (2^n_refinements);
            cells_y = ones(n_regions_y, 1) * (2^n_refinements);
            
            % 4. Instantiate the 2D Cartesian spatial mesh
            mesh = Mesh_2D_FDM(region_lengths_x, cells_x, region_lengths_y, cells_y);
            
            % 5. Vectorially upscale the 17x17 region map to a cell-by-cell material matrix
            % This step matches the resolution requirements of Solver_2D_FDM
            cell_materials = repelem(region_materials, cells_y, cells_x);
            
            % 6. Instantiate the 2D material library using the upscaled matrix
            nuclear_data = NuclearData_2D_FDM(cell_materials, D_lib, sigma_a_lib, ...
                                              nu_sigf_lib, chi_lib, sigma_s_lib);
                                          
            % 7. Instantiate the 2D system solver
            solver = Solver_2D_FDM(mesh, nuclear_data, bc);
        otherwise
            error("Unsupported dimension = %d. FDM solver supports 1D or 2D configurations.", dimension);
    end
    % System assembly and numerical solution
    fprintf("Assembling global matrices (A and B operators)...\n");
    solver = solver.assembleMatrices();
    
    fprintf("Compute the solution for the generalized eigenvalue problem %d mode(s)...\n", n_eigenvalues);
    solver = solver.solveEigenvalues(n_eigenvalues);
    
    % Display k-eff eigenvalues inside the terminal window
    solver.displayProblem();
    % Establish absolute paths to your explicit subfolders relative to root_path
    var_dir = fullfile(root_path, 'results', 'variables');
    fig_dir = fullfile(root_path, 'results', 'figures');
    % Defensive Architecture: Automatically generate folders if they don't exist yet
    if ~isfolder(var_dir), mkdir(var_dir); end
    if ~isfolder(fig_dir), mkdir(fig_dir); end
    % Extract the clean filename (e.g., "result") discarding any old extensions
    [~, out_name, ~] = fileparts(output_filename);
    if isempty(out_name)
        out_name = "fdm_result"; % Safe fallback name if config parameter is empty
    end
% --- Handle Graphical Export (results/figures) ---
    if output_flag
        fprintf("Rendering spatial neutron flux distributions...\n");
        
        % Loop through each requested eigenvalue/mode
        for m = 1:n_eigenvalues
            % Clear any previous figures to avoid cross-contamination
            close all;
            
            % Call the solver plotting method for the current eigenmode
            solver.plotPhi(m);
            
            % Automatically find all figure windows generated by this mode call
            fig_handles = findobj('Type', 'figure');
            
            % Loop through each opened window to save it uniquely
            for f = 1:length(fig_handles)
                h_curr = fig_handles(f);
                
                % Extract the window title (e.g., "Flux — Group 1, Mode 2")
                fig_title = get(h_curr, 'Name');
                
                % Sanitize the window title string to be a safe filename
                clean_title = regexprep(fig_title, '[^\w]', '_'); 
                clean_title = regexprep(clean_title, '_+', '_'); % Remove double underscores
                
                % Fallback naming convention if the window has no title string
                if isempty(clean_title)
                    clean_title = sprintf('mode_%d_plot_%d', m, f);
                end
                
                % Construct an absolute path with a completely unique identifier
                fig_export_path = char(fullfile(fig_dir, out_name + "_" + clean_title + ".png"));
                
                % Force window focus and apply generic presentation layouts
                figure(h_curr);
                grid on;
                
                % Save the image to your figures directory
                saveas(h_curr, fig_export_path);
                fprintf("Graphical visualization successfully saved to: '%s'\n", fig_export_path);
            end
        end
        % Close windows to clean up your MATLAB desktop layout workspace
        close all; 
    else
        fprintf("Plotting skipped (output_flag = false).\n");
    end
    
    % --- Handle Variable Export (results/variables) ---
    % Construct absolute script filename mapping inside the dedicated variables folder
    script_filename = char(fullfile(var_dir, out_name + ".m"));
    
    % Open text file stream for writing ('w')
    fid = fopen(script_filename, 'w');
    if fid == -1
        error("Core Engine Error: Unable to build variables file path at '%s'. Check permissions.", script_filename);
    end
    
    % Extract arrays from object attributes
    keff_vals = solver.keff;
    phi_vals  = solver.phi;
    
    % Track structural array bounds explicitly
    num_spatial_cells = size(phi_vals, 1);
    num_energy_groups = size(phi_vals, 2);
    num_eigenmodes    = size(phi_vals, 3);
    % Write Clear Descriptive Header
    fprintf(fid, '%% =========================================================================\n');
    fprintf(fid, '%% EXPORTED REACTOR SIMULATION WORKSPACE VARIABLES\n');
    fprintf(fid, '%% Generated on: %s\n', datestr(now));
    fprintf(fid, '%% Configuration Target Path: results/variables/%s.m\n', out_name);
    fprintf(fid, '%% =========================================================================\n\n');
    % Write Eigenvalues Vector (keff)
    fprintf(fid, '%% --- CRITICALITY EIGENVALUES [num_modes x 1] ---\n');
    fprintf(fid, 'keff = [\n');
    for m = 1:length(keff_vals)
        fprintf(fid, '    %.10f;\n', keff_vals(m));
    end
    fprintf(fid, '];\n\n');
    % Write Multi-Group Flux Matrix Slices (phi)
    fprintf(fid, '%% --- SPATIAL NEUTRON FLUX PROFILE ARRAY [Cells x Groups x Modes] ---\n');
    fprintf(fid, 'phi = zeros(%d, %d, %d);\n\n', num_spatial_cells, num_energy_groups, num_eigenmodes);
    for m = 1:num_eigenmodes
        fprintf(fid, '%% --- EIGENMODE OVERLAY: MODE %d ---\n', m);
        for g = 1:num_energy_groups
            fprintf(fid, 'phi(:, %d, %d) = [\n', g, m);
            for c = 1:num_spatial_cells
                fprintf(fid, '    %.8e;\n', phi_vals(c, g, m));
            end
            fprintf(fid, '];\n\n');
        end
    end
    % Safe termination of low-level data streams
    fclose(fid);
    fprintf("Physical variables successfully serialized to:  '%s'\n", script_filename); 
else
    error("Invalid method '%s'. Must be 'FEM' or 'FDM'.", method);
end

% =========================================================================
% AUXILIAR FUNCTIONS FEM
% =========================================================================

function params = build_params(energy_groups, n_refinements, fe_degree, perturbation_function, ...
                               boundary_conditions, output_filename, output_flag, ...
                               transient, time_delta, time_end, ...
                               region_lengths, cells_per_region, region_materials)
    params = struct();
    params.G = energy_groups;
    params.grado_l = fe_degree;
    params.n_refinements = n_refinements;
    params.tipo_pert = perturbation_function;
    params.boundary_conditions = boundary_conditions;
    params.output_filename = fullfile("results", "variables", output_filename);
    params.output_flag = output_flag;
    params.transient = transient;
    if transient
        params.time_delta = time_delta;
        params.t_total = time_end;
    end

    % Transform region data to FEM cells
    params.N = sum(cells_per_region);
    params.tam_celdas = zeros(1, params.N);
    params.materiales = zeros(1, params.N);
    idx = 1;
    for r = 1:length(region_lengths)
        cell_size = region_lengths(r) / cells_per_region(r);
        for c = 1:cells_per_region(r)
            params.tam_celdas(idx) = cell_size;
            params.materiales(idx) = region_materials(r);
            idx = idx + 1;
        end
    end
    params.L = sum(params.tam_celdas);
end

function params = merge_bank_data(params, bank_data)
    % Copy physical fields from bank (FDM -> FEM indexation)
    params.D_0 = bank_data.D_lib.';
    params.sigma_a_0 = bank_data.sigma_a_lib.';
    params.nu_sigma_f_0 = bank_data.nu_sigf_lib.';
    params.sigma_s_0 = permute(bank_data.sigma_s_lib, [2,3,1]);
    if size(bank_data.chi_lib, 1) > 1
        params.chi = bank_data.chi_lib(1,:).';
    else
        params.chi = bank_data.chi_lib.';
    end
    params.beta = bank_data.beta_prec;
    params.lambda = bank_data.lambda_prec;
    params.v = bank_data.v;

    % Compatibility check (number of materials)
    n_mat_bank = size(bank_data.D_lib, 1);
    if max(params.materiales) > n_mat_bank
        error('Incompatibility: materiales contains index %d, but bank has only %d material(s).', ...
              max(params.materiales), n_mat_bank);
    end
end

function params = time_functions(params)
    % Build default constant time functions
    params.D_func = cell(params.G, size(params.D_0, 2));
    params.sigma_a_func = cell(params.G, size(params.D_0, 2));
    params.nu_sigma_f_func = cell(params.G, size(params.D_0, 2));
    params.sigma_s_func = cell(params.G, params.G, size(params.D_0, 2));

    for g = 1:params.G
        for mat = 1:size(params.D_0, 2)
            params.D_func{g,mat} = @(t) params.D_0(g,mat);
            params.sigma_a_func{g,mat} = @(t) params.sigma_a_0(g,mat);
            params.nu_sigma_f_func{g,mat} = @(t) params.nu_sigma_f_0(g,mat);
            for h = 1:params.G
                if params.sigma_s_0(h,g,mat) ~= 0
                    params.sigma_s_func{h,g,mat} = @(t) params.sigma_s_0(h,g,mat);
                end
            end
        end
    end
end

function params = apply_perturbation(params, transient, perturbation_function, ...
                                     material_changing, slope)
    if transient && ~strcmp(perturbation_function, 'none')
        if isscalar(slope)
            slope = repmat(slope, size(material_changing));
        elseif length(slope) ~= length(material_changing)
            error('Slope must be scalar or have same length as material_changing');
        end

        switch lower(perturbation_function)
            case "linear"
                time_func = @(t, sl) sl * t;
            case "step"
                time_func = @(t, sl) sl * (t > 0);
            otherwise
                error('Type of perturbation "%s" not recognised', perturbation_function);
        end

        for i = 1:length(material_changing)
            mat = material_changing(i);
            sl = slope(i);
            for g = 1:params.G
                params.sigma_a_func{g,mat} = @(t) params.sigma_a_0(g,mat) + time_func(t, sl);
                params.nu_sigma_f_func{g,mat} = @(t) params.nu_sigma_f_0(g,mat) + time_func(t, sl);
            end
        end
    end
end

function [transitorio, t, y] = run_fem_solver(params, output_flag, Solver)
    addpath(fullfile(fileparts(mfilename('fullpath')), "FEM"));
    switch Solver
        case "fem_2gTransPC"
            [transitorio, t, y] = fem_2gTransPC(params, output_flag);
        case "fem_gTransPC"
            [transitorio, t, y] = fem_gTransPC(params, output_flag);
        case "fem_gTrans"
            [transitorio, t, y] = fem_gTrans(params, output_flag);
        otherwise
            error('Solver FEM "%s" not recognised', Solver);
    end
end

function [transitorio, t, y] = run_stationary(params, output_flag)
    malla = Malla1D(params.L, params.N, params.grado_l, params.tam_celdas);
    materiales_obj = Materiales1Dg(params.G, params.materiales, ...
        params.D_0, params.sigma_a_0, params.nu_sigma_f_0, ...
        params.sigma_s_0, params.chi);
    elemento = ElementoFinito(params.grado_l);
    estacionario = ProblemaDifusion1Dg(malla, materiales_obj, elemento);
    estacionario = estacionario.ensamblar_matrices_fem();
    estacionario = estacionario.aplicar_cc();
    estacionario = estacionario.resolver_autovalor(1);
    transitorio = estacionario;
    t = 0;
    y = estacionario.phi_inc(:, 1);

    if output_flag
        figure('Name', ['Estacionario: ', params.nombre]);
        x_nodos = malla.x_nodos;
        G = params.G;
        n_dof = length(x_nodos);
        hold on;
        colores = jet(G);
        for g = 1:G
            idx = (g-1)*n_dof + (1:n_dof);
            plot(x_nodos, y(idx), 'Color', colores(g,:), 'LineWidth', 2, ...
                 'DisplayName', ['Grupo ', num2str(g)]);
        end
        xlabel('x (cm)');
        ylabel('Flujo');
        legend('Location', 'best');
        grid on;
        title(['Flujo estacionario - ', params.nombre]);
        saveas(gcf, fullfile('results/figures', ['Estacionario_', params.nombre, '.eps']), 'epsc');
    end
end