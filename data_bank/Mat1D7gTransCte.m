% Mat1D7gTransCte
G = 7;
D_lib = [1.5, 1.2, 0.9, 0.6, 0.4, 0.3, 0.2];
sigma_a_lib = [0.0025, 0.003, 0.004, 0.006, 0.01, 0.02, 0.04];
nu_sigf_lib = [0.002, 0.0015, 0.001, 0.0005, 0.0002, 0.0001, 0];
sigma_s_lib = zeros(1, 7, 7);
sigma_s_lib(1, 2, 1) = 0.025;
sigma_s_lib(1, 3, 2) = 0.020;
sigma_s_lib(1, 4, 3) = 0.015;
sigma_s_lib(1, 5, 4) = 0.010;
sigma_s_lib(1, 6, 5) = 0.008;
sigma_s_lib(1, 7, 6) = 0.005;
chi_lib = [0.6, 0.3, 0.1, 0, 0, 0, 0];
beta_prec = [];
lambda_prec = [];
v = [3e7, 2e7, 1e7, 5e6, 2e6, 5e5, 2e5];