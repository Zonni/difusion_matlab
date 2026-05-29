function mat_obj = Materiales1Dg(G, materiales_vec, D_0, sigma_a_0, nu_sigma_f_0, sigma_s_0, chi)
    % Crea un objeto de materiales para G grupos de energia
    % D_0, sigma_a_0, nu_sigma_f_0: matriz (G x n_mat)
    % sigma_s_0: array (G x G x n_mat)
    % chi: vector (G x 1)
    % materiales_vec: vector con el indice del material en cada celda (1 x N)
    mat_obj.G = G;
    mat_obj.D = D_0;
    mat_obj.sigma_a = sigma_a_0;
    mat_obj.nu_sigma_f = nu_sigma_f_0;
    mat_obj.sigma_s = sigma_s_0;
    mat_obj.chi = chi;
    mat_obj.ind_materiales = materiales_vec;
end