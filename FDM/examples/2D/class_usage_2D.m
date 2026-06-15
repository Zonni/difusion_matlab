% EXAMPLE_2D_MULTIGROUP
% Usage example for the 2D multigroup neutron diffusion solver.
%
% Geometry (7 x 7 cells, 10 cm each = 70 x 70 cm total):
%
%   0  1  1  1  1  1  0
%   1  1  2  2  2  1  1
%   1  2  2  2  2  2  1
%   1  2  2  2  2  2  1
%   1  2  2  2  2  2  1
%   1  1  2  2  2  1  1
%   0  1  1  1  1  1  0
%
%   0 = outside reactor (vacuum boundary applied to adjacent active faces)
%   1 = reflector  (no fission)
%   2 = fuel       (fissions in group 2 / thermal)
%
% Two energy groups:
%   Group 1 — fast neutrons
%   Group 2 — thermal neutrons
%
% Fission neutrons are born in group 1 only (chi = [1, 0]).
% Downscattering from group 1 to group 2 is the only scattering channel.

clearvars; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'src', '2D'));

%% ========================================================================
%  1. MESH
%  Each cell is 10 cm wide and 10 cm tall.  The full domain is 70 x 70 cm.
% =========================================================================

cell_size = 10;           % cm
n_cells   = 7;

delta_x = cell_size * ones(n_cells, 1);   % (7 x 1)
delta_y = cell_size * ones(n_cells, 1);   % (7 x 1)

mesh = Mesh_2D_FDM(delta_x, delta_y);
mesh.displayMesh();
mesh.plotMesh();

%% ========================================================================
%  2. NUCLEAR DATA
%  Material library (2 materials x 2 groups).
%
%  Cross-section layout: rows = materials, cols = energy groups.
%  sigma_s_lib dimensions: (material, g_from, g_to).
%  Only sigma_s(mat, 1, 2) != 0 — downscattering fast->thermal.
% =========================================================================

% --- Diffusion coefficients (cm) ----------------------------------------
%               group 1    group 2
D_lib      = [ 1.446,     0.276;  ...   % material 1: reflector
               1.268,     0.367 ];      % material 2: fuel

% --- Absorption cross sections (1/cm) -----------------------------------
sigma_a_lib = [ 0.0002,   0.0100;  ...  % reflector
                0.0082,   0.0730 ];     % fuel

% --- Nu * fission cross sections (1/cm) ---------------------------------
%  Reflector has no fission; fuel fissions only in the thermal group.
nu_sigf_lib = [ 0,        0;       ...  % reflector
                0,        0.1660 ];     % fuel

% --- Fission spectrum ---------------------------------------------------
%  All fission neutrons born as fast (group 1).
chi_lib     = [ 1,        0;       ...  % reflector (unused — nu_sigf = 0)
                1,        0 ];          % fuel

% --- Scattering matrix (material, g_from, g_to) -------------------------
%  Size must be (n_materials x n_groups x n_groups) = (2 x 2 x 2).
%  Only downscattering (1->2) is non-zero.
sigma_s_lib             = zeros(2, 2, 2);
sigma_s_lib(1, 1, 2)    = 0.0200;   % reflector: fast->thermal
sigma_s_lib(2, 1, 2)    = 0.0170;   % fuel:      fast->thermal

% --- Material grid (Ny x Nx = 7 x 7) -----------------------------------
%  Row index = y direction (top to bottom), col index = x direction.
region_materials = [ ...
    0,  1,  1,  1,  1,  1,  0 ; ...
    1,  1,  2,  2,  2,  1,  1 ; ...
    1,  2,  2,  2,  2,  2,  1 ; ...
    1,  2,  2,  2,  2,  2,  1 ; ...
    1,  2,  2,  2,  2,  2,  1 ; ...
    1,  1,  2,  2,  2,  1,  1 ; ...
    0,  1,  1,  1,  1,  1,  0 ];

nuclear_data = NuclearData_2D(region_materials, D_lib, sigma_a_lib, ...
                              nu_sigf_lib, chi_lib, sigma_s_lib);
nuclear_data.displayMaterials();

%% ========================================================================
%  3. SOLVER
% =========================================================================

solver = Solver_2D_FDM(mesh, nuclear_data);
solver = solver.assembleMatrices();
solver = solver.solveEigenvalues(1);   % compute fundamental mode only

solver.displayProblem();
solver.plotPhi(1);                     % plot fundamental flux, both groups