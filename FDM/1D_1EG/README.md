# 1D 1-Group Neutron Diffusion Solver

This project implements a numerical solver in MATLAB to analyze steady-state reactor physics problems. It calculates the fundamental multiplication factor ($k_{eff}$) and the neutron flux profile ($\Phi(x)$) for heterogeneous slab geometries.

Currently, the solver uses the **Mesh-Centered Finite Difference Method (FDM)**, roughly following the mathematical formulations established in Alain Hébert's _Applied Reactor Physics_.

## Key Features

-   **Object-Oriented Architecture:** Clean separation of concerns between Geometry (`Mesh_1D_FDM`), Physics (`Materials`), and the Mathematical formulation (`Solver_1D_1EG_FDM`).
    
-   **Heterogeneous Domain Support:** Easily define multi-region reactors (e.g., a fissile core surrounded by non-multiplying reflectors) with distinct cross-sections.
    
-   **Eigenvalue Solver:** Computes $k_{eff}$ and both fundamental and harmonic flux modes using MATLAB's sparse matrix generalized eigenvalue solver.
    

## Mathematical Formulation

The solver discretizes the 1D 1-group steady-state neutron diffusion equation:

$$-\frac{d}{dx} \left( D(x) \frac{d\phi(x)}{dx} \right) + \Sigma_a(x) \phi(x) = \frac{1}{k_{eff}} \nu\Sigma_f(x) \phi(x)$$

Using a finite-volume integration approach over discrete cells, this is cast into the algebraic eigenvalue problem:

$$\mathbb{A} \vec{\Phi} = \frac{1}{k_{eff}} \mathbb{F} \vec{\Phi}$$

Where $\mathbb{A}$ is the strictly positive-definite removal/transport matrix (accounting for absorption and leakage) and $\mathbb{F}$ is the positive fission production matrix.

## File Overview

Currently, the project is housed in a flat directory structure for rapid development and testing. The core classes include:

-   `Mesh_1D_FDM.m`: Defines the spatial domain, discretizing physical regions into computational cells.
    
-   `clean_materials.m`: Acts as the library for macroscopic cross-sections ($\Sigma_a, \nu\Sigma_f$) and diffusion coefficients ($D$).
    
-   `Solver_1D_1EG_FDM.m`: Assembles the linear system matrices and solves for the eigenvalues and flux eigenvectors.

