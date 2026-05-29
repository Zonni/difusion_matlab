
# Mesh_1D_FDM - Implementation Notes

## Mesh_1D_FDM: domain discretization vectorization
```MATLAB
x_left = obj.region_boundaries(obj.cell_region_idx);
temp_before = repelem([0; cumsum(obj.cells_per_region(1:end-1))], obj.cells_per_region);
cells_before = temp_before(:);

j_local = (1:obj.num_cells)' - cells_before;
obj.cell_centers = x_left + obj.delta_x .* (j_local - 0.5);
```

For every region $i$, we want to obtain a uniform mesh of centered points. This is,

$$x_j = region\_boundaries(i-1)+ \Delta x_i \ (j - 1/2), \\ \forall j \in \{1, \dots, cells\_per\_region(i)\}.$$

In order to vectorize this procedure, let us define an index $k = 1, \dots, num\_cells$, that maps each number to the corresponding $k$-th node in the mesh $x_k$. Now, given a region $i$ and a node in this subregion of the mesh $x_j$, we have the following relation: 

$$j = k - cells\_before(i),$$

where $cells\_before(i)$ is the total number of cells in the regions prior to region $i$. We can think of it as "truncating" $k$ in region $i$ so to make it start at 1.

This justifyes the vectorial representation

$$cell\_centers = x\_left+ delta\_x \ (j\_local - 1/2).$$
