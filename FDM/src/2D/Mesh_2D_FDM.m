classdef Mesh_2D_FDM
    %MESH_2D_FDM Geometry manager for 2D Cartesian discretization.
    %   Creates and manages the spatial mesh for a 2D finite difference
    %   method solver. Regions are defined independently along X and Y;
    %   each region is subdivided into a uniform number of cells, exactly
    %   as in Mesh_1D_FDM.

    properties
        % X-Axis
        region_lengths_x    (:,1) double    % Length of each X region (cm)
        cells_per_region_x  (:,1) double    % Number of cells per X region
        num_cells_x         (1,1) double    % Total number of cell columns
        cell_centers_x      (:,1) double    % X coordinate of cell centres (Nx x 1)
        cell_region_idx_x   (:,1) double    % X region index each column belongs to
        delta_x             (:,1) double    % Width of each cell column (Nx x 1)
        region_boundaries_x (:,1) double    % Cumulative X boundaries (Nx_reg+1 x 1)

        % Y-Axis
        region_lengths_y    (:,1) double    % Length of each Y region (cm)
        cells_per_region_y  (:,1) double    % Number of cells per Y region
        num_cells_y         (1,1) double    % Total number of cell rows
        cell_centers_y      (:,1) double    % Y coordinate of cell centres (Ny x 1)
        cell_region_idx_y   (:,1) double    % Y region index each row belongs to
        delta_y             (:,1) double    % Height of each cell row (Ny x 1)
        region_boundaries_y (:,1) double    % Cumulative Y boundaries (Ny_reg+1 x 1)

        % Global
        num_cells_total     (1,1) double
    end

    methods
        function obj = Mesh_2D_FDM(region_lengths_x, cells_per_region_x, ...
                                   region_lengths_y, cells_per_region_y)
            %MESH_2D_FDM Construct an instance of this class.
            %
            %   Inputs:
            %       region_lengths_x   - (Nx_reg x 1) Length of each X region (cm).
            %       cells_per_region_x - (Nx_reg x 1) Number of uniform cells per X region.
            %       region_lengths_y   - (Ny_reg x 1) Length of each Y region (cm).
            %       cells_per_region_y - (Ny_reg x 1) Number of uniform cells per Y region.
            %
            %   Example — 3-region X axis, 2-region Y axis:
            %       mesh = Mesh_2D_FDM([10;20;10], [2;4;2], [15;30], [3;6]);

            % Security size checks
            if ~isequal(size(region_lengths_x), size(cells_per_region_x))
                error("region_lengths_x and cells_per_region_x must have the same size.");
            end
            if ~isequal(size(region_lengths_y), size(cells_per_region_y))
                error("region_lengths_y and cells_per_region_y must have the same size.");
            end
            if any(region_lengths_x(:) <= 0) || any(region_lengths_y(:) <= 0)
                error("All region lengths must be strictly positive.");
            end
            if any(cells_per_region_x(:) < 1) || any(cells_per_region_y(:) < 1)
                error("Each region must contain at least one cell.");
            end

            % --- X axis --------------------------------------------------
            obj.region_lengths_x    = region_lengths_x(:);
            obj.cells_per_region_x  = cells_per_region_x(:);
            obj.region_boundaries_x = [0; cumsum(obj.region_lengths_x)];
            obj.num_cells_x         = sum(obj.cells_per_region_x);

            obj.cell_region_idx_x = repelem( ...
                (1:length(obj.region_lengths_x))', obj.cells_per_region_x);

            obj.delta_x = obj.region_lengths_x(obj.cell_region_idx_x) ./ ...
                          obj.cells_per_region_x(obj.cell_region_idx_x);

            x_left         = obj.region_boundaries_x(obj.cell_region_idx_x);
            cells_before_x = repelem( ...
                [0; cumsum(obj.cells_per_region_x(1:end-1))], obj.cells_per_region_x);
            j_local_x      = (1:obj.num_cells_x)' - cells_before_x(:);
            obj.cell_centers_x = x_left + obj.delta_x .* (j_local_x - 0.5);

            % --- Y axis --------------------------------------------------
            obj.region_lengths_y    = region_lengths_y(:);
            obj.cells_per_region_y  = cells_per_region_y(:);
            obj.region_boundaries_y = [0; cumsum(obj.region_lengths_y)];
            obj.num_cells_y         = sum(obj.cells_per_region_y);

            obj.cell_region_idx_y = repelem( ...
                (1:length(obj.region_lengths_y))', obj.cells_per_region_y);

            obj.delta_y = obj.region_lengths_y(obj.cell_region_idx_y) ./ ...
                          obj.cells_per_region_y(obj.cell_region_idx_y);

            y_left         = obj.region_boundaries_y(obj.cell_region_idx_y);
            cells_before_y = repelem( ...
                [0; cumsum(obj.cells_per_region_y(1:end-1))], obj.cells_per_region_y);
            j_local_y      = (1:obj.num_cells_y)' - cells_before_y(:);
            obj.cell_centers_y = y_left + obj.delta_y .* (j_local_y - 0.5);

            obj.num_cells_total = obj.num_cells_x * obj.num_cells_y;
        end

        % -----------------------------------------------------------------
        function displayMesh(obj)
            %DISPLAYMESH Prints a summary of the 2D mesh geometry.
            fprintf('=================================\n        2D MESH DATA        \n=================================\n');
            fprintf('Total Cells : %d  (Nx = %d, Ny = %d)\n', ...
                    obj.num_cells_total, obj.num_cells_x, obj.num_cells_y);
            fprintf('Total Length: X = %.2f cm,  Y = %.2f cm\n\n', ...
                    obj.region_boundaries_x(end), obj.region_boundaries_y(end));

            fprintf('--- X Regions ---\n');
            Region_X = (1:length(obj.region_lengths_x))';
            Length_X = obj.region_lengths_x;
            Cells_X  = obj.cells_per_region_x;
            disp(table(Region_X, Length_X, Cells_X));

            fprintf('--- Y Regions ---\n');
            Region_Y = (1:length(obj.region_lengths_y))';
            Length_Y = obj.region_lengths_y;
            Cells_Y  = obj.cells_per_region_y;
            disp(table(Region_Y, Length_Y, Cells_Y));
        end

        % -----------------------------------------------------------------
        function plotMesh(obj)
            %PLOTMESH Generates a 2D plot of the finite difference grid,
            %   showing internal cell boundaries (light grey) and physical
            %   region boundaries (thick black).

            figure('Name', '2D FDM Mesh Grid', 'Color', 'w');
            hold on; axis equal; box on;
            title('2D FDM Spatial Mesh');
            xlabel('X Position (cm)');
            ylabel('Y Position (cm)');

            % Internal cell boundaries (light grey)
            x_faces = [0; cumsum(obj.delta_x)];
            y_faces = [0; cumsum(obj.delta_y)];

            for i = 1:numel(x_faces)
                plot([x_faces(i), x_faces(i)], ...
                     [0, obj.region_boundaries_y(end)], ...
                     'Color', [0.85 0.85 0.85]);
            end
            for j = 1:numel(y_faces)
                plot([0, obj.region_boundaries_x(end)], ...
                     [y_faces(j), y_faces(j)], ...
                     'Color', [0.85 0.85 0.85]);
            end

            % Physical region boundaries (thick black)
            for i = 1:numel(obj.region_boundaries_x)
                xb = obj.region_boundaries_x(i);
                plot([xb, xb], [0, obj.region_boundaries_y(end)], ...
                     'k', 'LineWidth', 1.8);
            end
            for j = 1:numel(obj.region_boundaries_y)
                yb = obj.region_boundaries_y(j);
                plot([0, obj.region_boundaries_x(end)], [yb, yb], ...
                     'k', 'LineWidth', 1.8);
            end

            xlim([0, obj.region_boundaries_x(end)]);
            ylim([0, obj.region_boundaries_y(end)]);
            hold off;
        end
    end
end