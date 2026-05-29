classdef Mesh_1D_FDM
    %MESH_1D_FDM Geometry manager for 1D slab discretization.
    %   Creates and manages the spatial mesh for a 1D finite difference 
    %   method solver, given region lengths and discretizations.
    properties
        region_lengths    (:,1) double
        cells_per_region  (:,1) double
        num_cells         (1,1) double
        cell_centers      (:,1) double
        cell_region_idx   (:,1) double
        delta_x           (:,1) double
        region_boundaries (:,1) double
    end

    methods
        function obj = Mesh_1D_FDM(region_lengths, cells_per_region)
            %MESH_1D_FDM Construct an instance of this class.
            %
            %   Inputs:
            %       region_lengths   - (N x 1) Array of lengths for each physical region (cm).
            %       cells_per_region - (N x 1) Array of the number of uniform cells per region.
            %
            %   Outputs:
            %       obj              - Initialized Mesh_1D_FDM object.
            
            % Security size check
            if ~isequal(size(region_lengths), size(cells_per_region))
                error("region_lengths and cells_per_region must have the same size.")
            end

            obj.region_lengths = region_lengths(:);
            obj.cells_per_region = cells_per_region(:);
            obj.region_boundaries = [0; cumsum(obj.region_lengths)];
            obj.num_cells = sum(obj.cells_per_region);
            
            % Create a uniform mesh in each region
            cell_region_idx = repelem((1:length(region_lengths))', cells_per_region);
            obj.cell_region_idx = cell_region_idx(:);
            obj.delta_x = obj.region_lengths(obj.cell_region_idx) ./ obj.cells_per_region(obj.cell_region_idx);
            
            x_left = obj.region_boundaries(obj.cell_region_idx);
            temp_before = repelem([0; cumsum(obj.cells_per_region(1:end-1))], obj.cells_per_region);
            cells_before = temp_before(:); 
            
            j_local = (1:obj.num_cells)' - cells_before;
            obj.cell_centers = x_left + obj.delta_x .* (j_local - 0.5);        
        end
        
        function displayMesh(obj)
            %DISPLAYMESH Prints a summary of the 1D mesh geometry.
            %
            %   Inputs:
            %       obj - Mesh_1D_FDM object.
            %
            %   Outputs:
            %       None.
            fprintf('======================\n    MESH DATA    \n======================\n\n');
            fprintf('Total Length: %.2f cm\n', obj.region_boundaries(end));
            fprintf('Total Cells:  %d\n',      obj.num_cells);
            fprintf('Regions:      %d\n\n',    length(obj.region_lengths));
            Region = (1:length(obj.region_lengths))';
            Length = obj.region_lengths;
            Cells  = obj.cells_per_region;
            disp(table(Region, Length, Cells));
        end
    end
end