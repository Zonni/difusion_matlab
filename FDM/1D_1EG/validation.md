# Code Validation

**Numerical Methods**
* Jaime's FDM
* Jaime's FEM
* Ruben's FDM

## Homogeneous reactor
**Material Properties**
* The whole domain is Fuel (index 2).

In the script, we start by computing the analytical $k_{eff}$. For a bare, homogeneous 1D slab of length $L$, the analytical effective multiplication factor for the fundamental mode ($n=1$) can be obtained as:

$$k_{eff} = \frac{\nu\Sigma_f}{\Sigma_a + D\left(\frac{\pi}{L}\right)^2}.$$

Then, we define the analytical solution of the flux as an anonymous function using the appropiate formula:

$$sin(\frac{\pi x}{L}),$$ 

satisfying the Dirichlet boundary conditions ($\phi(0) = \phi(L) = 0$).

Instead of looping over mesh sizes per se, we loop over the "desired degrees of freedom" (`target_dof`) of the systems, because Jaime's mesh implementation takes a given number of nodes in the mesh and adds more depending on the degree of the FEM. In order to relate the `target_dof` with the number of nodes in Jaime's mesh, we use the formula:

$$N_{FEM} = \max\left(1, \text{round}\left(\frac{\text{target\_dof} - 1}{\text{FEM\_DEGREE}}\right)\right),$$

based on his implementation.

For all methods, we compute their respective numerical approximation. Then, we use the method's spatial mesh array to evaluate the analytical solution and normalize said analytical solution vector (because the numerical vector fluxes are normalized in every implementation). Then, using MATLAB's `rmse` function, we get the global flux error. The Root Mean Square Error (RMSE) quantifies the deviation of the numerical profile $\phi_{num}$ from the analytical profile $\phi_{ana}$ across all nodes:

$$RMSE = \sqrt{\frac{1}{N_{nodes}}\sum_{i=1}^{N_{nodes}}\left(\phi_{ana}(x_i) - \phi_{num}(x_i)\right)^2}.$$

It is not exactly the formula recquired by the research team, but it is the one Jaime used and the one with which I could verify his results.

In the end, we compute the $k_{eff}$ error using this formula, which expresses the absolute difference in **pcm (por cien mil)**:
$$\epsilon_k = |k_{eff}^{ana} - k_{eff}^{num}| \times 10^5,$$

and print on screen tables and produce plots.

---

## Heterogeneous reactor (Reflector-Core-Reflector)
**Material Properties**
* The domain is divided into three regions: Reflector (index 1), Fuel (index 2), and Reflector (index 1).
*The reflector regions have no fission production ($\nu\Sigma_f = 0$).

Unlike the homogeneous case, the heterogeneous analytical solution cannot be expressed as a simple sine wave due to the varying material properties and the physical discontinuities at the interfaces. 

To compute the exact analytical solution, we utilize a symmetric coordinate transformation. By shifting the origin ($x=0$) to the exact center of the core, we force the flux profile to be perfectly symmetric, allowing us to define the flux in the core using a simple cosine function and the flux in the reflector using a hyperbolic sine function that decays to zero at the outer boundary.



**1. Initial Guess via Extrapolation Distance**
To solve for $k_{eff}$, the root-finding algorithm (`fzero`) requires an initial guess. We approximate the reflected core as a slightly larger "bare" core by adding the extrapolation distance ($d_{ext}$) to the physical core width. 



$$d_{ext} \approx 2.13 D_c$$

We then calculate the initial guess using the bare reactor formula:
$$k_{bare} = \frac{\nu\Sigma_f}{D_c B_g^2 + \Sigma_{ac}}$$
where $B_g = \frac{\pi}{2(a_{core} + d_{ext})}$ is the geometric buckling of the extrapolated bare core.

**2. The Criticality Equation**
[cite_start]We apply the physical interface conditions at the core-reflector boundary: continuity of flux and continuity of current (Fick's Law)[cite: 270]. Dividing the current equation by the flux equation eliminates the amplitude coefficients and yields the exact transcendental criticality equation:
$$D_c B_c \tan(B_c a_{core}) = D_r \kappa_r \coth(\kappa_r b_{refl})$$
Where:
* $B_c = \sqrt{\frac{\frac{\nu\Sigma_f}{k_{eff}} - \Sigma_{ac}}{D_c}}$ is the material buckling of the core.
* [cite_start]$\kappa_r = \sqrt{\frac{\Sigma_{ar}}{D_r}}$ is the inverse diffusion length of the reflector[cite: 256].

MATLAB's `fzero` solves this equation to find the exact analytical $k_{eff}$.

**3. The Flux Profile**
Once $k_{eff}$ is known, the material buckling $B_c$ is fixed. Setting the arbitrary center maximum flux to $A = 1.0$, we determine the reflector amplitude coefficient $C$ via the flux continuity equation:
$$C = \frac{A \cos(B_c a_{core})}{\sinh(\kappa_r b_{refl})}$$

The exact analytical flux evaluated at any point $x$ (relative to the core center) is:
* **Core ($|x| \le a_{core}$):** $\phi(x) = A \cos(B_c x)$
* **Reflector ($|x| > a_{core}$):** $\phi(x) = C \sinh(\kappa_r(a_{core} + b_{refl} - |x|))$

During validation, because the numerical vectors (`phi_num`) from the FDM and FEM methods are natively normalized by their respective solvers, we scale the analytical flux to match the maximum absolute value of the numerical flux before executing the RMSE calculation. [cite_start]The $k_{eff}$ error is calculated using the standard pcm formula[cite: 308].
