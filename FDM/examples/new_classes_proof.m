clearvars;
close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));

region_lenghts_x = [0.5; 0.5; 0.5]

region_lenghts_y = [0.5; 0.5; 0.5]

mesh2d = Mesh_2D_FDM(region_lenghts_x, region_lenghts_y)
mesh2d.region_boundaries_x
mesh2d.region_boundaries_y
mesh2d.num_cells
mesh2d.cell_centers
mesh2d.delta_x
mesh2d.delta_y