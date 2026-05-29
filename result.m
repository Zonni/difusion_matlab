% =========================================================================
% NEUTRON DIFFUSION SOLVER SIMULATION EXPORT
% Generated automatically on: 26-May-2026 19:14:04
% Method: FDM | Total Spatial DOF: 6
% =========================================================================

% --- CRITICALITY EIGENVALUES ---
% Effective multiplication factor vector [num_eigenmodes x 1]
keff = [
    0.8045428360;
];

% --- SPATIAL FLUX MULTI-ARRAY PROFILES ---
% Pre-allocating complete tensor: [Cells x Energy Groups x Eigenmodes]
phi = zeros(6, 1, 1);

% ==========================================
% Eigenmode Mode 1 Profile Array 
% ==========================================
phi(:, 1, 1) = [
    0.00000000e+00;
    2.16179393e-01;
    1.00000000e+00;
    1.00000000e+00;
    2.16179393e-01;
    0.00000000e+00;
];

