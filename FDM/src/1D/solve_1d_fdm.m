function [k_eff, phi, full_mesh_coordinates] = solve_1d_fdm(mesh, nuclear_data, nm)
%SOLVE_1D_FDM Wrapper function for Mesh_1D_FDM and NuclearData_1D
%   Given the class instances, solves the associated diffusion problem.
    arguments
        mesh Mesh_1D_FDM
        nuclear_data NuclearData_1D
        nm (1,1) {mustBePositive, mustBeInteger} = 1
    end
    
    solver = Solver_1D_FDM(mesh, nuclear_data);

    solver = solver.assembleMatrices().solveEigenvalues(nm);
    k_eff = solver.k_eff;
    phi = solver.phi;
    full_mesh_coordinates = solver.full_mesh_coordinates;
end
