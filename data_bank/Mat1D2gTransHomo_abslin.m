% Mat1D2gTransHomo_abslin
G = 2;
D_lib = [1.32, 0.2772];
sigma_a_lib = [0.0026562, 0.071596];
nu_sigf_lib = [0.0074527, 0.08391];
sigma_s_lib = zeros(1, 2, 2);
sigma_s_lib(1, 2, 1) = 0.023106;
chi_lib = [1.0, 0.0];
beta_prec = [];
lambda_prec = [];
v = [1.27e7; 2.5e5];