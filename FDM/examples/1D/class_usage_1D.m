% This script provides an example on how to introduce data into the
% 1D classes.

clearvars;
close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', '1D'));

%%

% 1D 1EG Problem with 2 materials

% MESH DATA
region_lengths = [10; 20; 10];
cells_per_region = [1; 2; 1];

% NUCLEAR DATA
region_materials = [1;2;1];
D_lib = [1.446; 0.776]; % in cm
sigma_a_lib = [0.0077; 0.0244]; % in 1/cm
nu_sigf_lib = [0; 0.0260]; % Mat 1 has no fission; 1/cm
chi_lib = [1; 1];
sigma_s_lib = zeros(2, 1, 1);

% Material grid and mesh initialization
mesh = Mesh_1D_FDM(region_lengths, cells_per_region);
mesh.displayMesh;

materials = NuclearData_1D(region_materials, D_lib, sigma_a_lib, nu_sigf_lib, chi_lib, sigma_s_lib);
materials.displayMaterials;

% Solver initialization
solver = Solver_1D_FDM(mesh, materials);

% Assembling the system and solving the eigenvalue problem
solver = solver.assembleMatrices.solveEigenvalues(1);
solver.displayProblem;
solver.plotPhi;