classdef mesh_2d_refinment
    %MESH_2D_FDM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        region_lengths_x
        region_lengths_y
        refinment_x
        refinment_y
        region_boundaries_x
        region_boundaries_y
        num_cells
        cell_centers
        cell_region_index
    end
    
    methods
        function obj = mesh_2d_refinment(region_lengths_x,region_lengths_y,refinment_x,refinment_y)
            %MESH_2D_FDM Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.region_lengths_x = region_lengths_x(:);
            obj.region_lengths_y = region_lengths_y(:);
            
            obj.refinment_x = refinment_x(:);
            obj.refinment_y = refinment_y(:);

            obj.region_boundaries_x = [0; cumsum(obj.region_lengths_x)];
            obj.region_boundaries_y = [0; cumsum(obj.region_lengths_y)];

            % Get the total number of cells:
            num_regions_x = length(obj.region_lengths_x);
            num_regions_y = length(obj.region_lengths_y);
            obj.num_cells = num_regions_x * num_regions_y;

            % Get the cell boundaries vectors
            cell_boundaries_x = zeros(length(obj.region_boundaries_x) + sum(refinment_x),1);
            for i = 1:length(obj.region_boundaries_x)
                if ~isequal(, size(cells_per_region))
            end
            
            % Get the coordinates of the cells
            obj.cell_centers = zeros(obj.num_cells, 2);

            for j = 1:num_regions_y
                obj.cell_centers(((j-1)*num_regions_x + 1):(j * num_regions_x),2) = obj.region_boundaries_y(j) + obj.region_lengths_y(j)/2
                for i = 1:num_regions_x
                    obj.region_boundaries_x(i) + obj.region_lengths_y(i)/2
                    obj.cell_centers((j-1)*num_regions_x + i,1) = obj.region_boundaries_x(i) + obj.region_lengths_y(i)/2
                end
            end

        end
    end
end

