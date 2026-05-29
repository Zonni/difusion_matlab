function [transitorio, t, y] = fem_gTrans(params, bool_graf)
    %% fem_gTrans
    % Resuelve un problema de difusion 1D transitorio con G grupos de energia
    % utilizando un integrador de ODE (ode15s)

    %% Estado estacionario
    malla = Malla1D(params.L, params.N, params.grado_l, params.tam_celdas);
    materiales_obj = Materiales1Dg(params.G, params.materiales, ...
        params.D_0, params.sigma_a_0, params.nu_sigma_f_0, ...
        params.sigma_s_0, params.chi);
    elemento = ElementoFinito(params.grado_l);

    estacionario = ProblemaDifusion1Dg(malla, materiales_obj, elemento);
    estacionario = estacionario.ensamblar_matrices_fem();
    estacionario = estacionario.aplicar_cc();
    estacionario = estacionario.resolver_autovalor(1);

    G = params.G;
    n_dof = malla.N * malla.grado_l + 1;
    Phi0 = estacionario.phi_inc(:, 1);
    keff0 = estacionario.keff(1);

    %% Transitorio
    n_mat = size(params.D_0, 2);
    % Forzar criticidad dividiendo las funciones de fision por keff0
    nu_sigma_f_n = cell(G, n_mat);
    for g = 1 : G
        for mat = 1 : n_mat
            nu_sigma_f_n{g, mat} = @(t) params.nu_sigma_f_func{g, mat}(t) / keff0;
        end
    end

    transitorio = ProblemaDifusion1DgTrans(malla, materiales_obj, elemento);
    transitorio.nombre_banco = params.nombre;
    transitorio.v = params.v(:);
    transitorio.beta = params.beta;
    transitorio.lambda = params.lambda;
    transitorio.n_precursor = length(params.beta);

    % Asignar funciones temporales (ya normalizadas)
    transitorio.D_func = params.D_func;
    transitorio.sigma_a_func = params.sigma_a_func;
    transitorio.nu_sigma_f_func = nu_sigma_f_n;
    transitorio.sigma_s_func = params.sigma_s_func;
    transitorio.keff0 = keff0;

    transitorio = transitorio.ensamblar_matrices_cte();
    y0 = transitorio.cond_ini(Phi0, 0);

    int_t = [0, params.t_total];
    [t, y] = transitorio.solver(int_t, y0, odeset());
    
    %% Resultados
    transitorio.Phi_t = y(:, 1 : G * n_dof);
    transitorio.Pot_num = transitorio.pot_norm(t, y);

    if bool_graf
        transitorio.graficar_t(t, y);
    end
end