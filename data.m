% --- GENERAL CONFIGURATION ---
% Numerical method (FEM or FDM)
method = "FEM";

% Number of eigenvalues required
n_eigenvalues = 1;

% Energy groups
energy_groups = 2;

% Geometrical parameters
mesh_size = 60;
%cell_pitch_x = 60* 5.00; A DATABANK
%cell_pitch_y = 40;

% --- INPUT OTIONS ---
bank = "Mat1D2gTransHomo_abslin";

% --- OUT OPTIONS ---
% Filename where will be written the output.
output_filename = "result.mat";

% True/false - Save graphs
output_flag = false;

% --- FDM RELATED CONFIGURATION ---

% Dimension of the problem (1 or 2)
dimension = 1;

% Number of global refinements
n_refinements = 0;

% Boundary Conditons (0 ZeroFlux) (1 Vaccum)
boundary_conditions = 0;

% --- FEM RELATED CONFIGURATION ---
%Materiales = [1 * ones(1, 14)]; A DATABANK

% Finite element degree		NOMES FEM
fe_degree = 1;

% TIME OPTIONS
transient = true;
Solver = "fem_gTransPC";

% Delta time of the iteration	NOMES SI FEM I TRANSIENT TRUE
time_delta = 0.001;
time_end = 10.0;

% Distribution of the instability: Flux_Distributed or Single_Material NOMES SI FEM I TRANSIENT
type_perturbation = "Ramp_Two_Mats";
material_changing = 2;

% Type of instability Constant | Ramp NOMES SI FEM I TRANSIENT
perturbation_function = "Ramp";
slope = -0.1;
