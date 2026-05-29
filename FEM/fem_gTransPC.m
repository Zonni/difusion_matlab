function [transitorio, t, y] = fem_gTransPC(params, bool_graf)
    %% fem_gTransPC
    % Resuelve el transitorio con el esquema semi‑implicito de ecuaciones (53) y (54)

    %% Paso 1: Estacionario
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

    %% Paso 2: Transitorio
    n_mat = size(params.D_0, 2);

    % Forzar criticidad
    nu_sigma_f_ren = cell(G, n_mat);
    for g = 1 : G
        for mat = 1:n_mat
            nu_sigma_f_ren{g, mat} = @(t) params.nu_sigma_f_func{g, mat}(t) / keff0;
        end
    end

    transitorio = ProblemaDifusion1DgTrans(malla, materiales_obj, elemento);
    transitorio.nombre_banco = params.nombre;
    transitorio.v = params.v(:);
    transitorio.beta = params.beta;
    transitorio.lambda = params.lambda;
    transitorio.n_precursor = length(params.beta);

    transitorio.D_func = params.D_func;
    transitorio.sigma_a_func = params.sigma_a_func;
    transitorio.nu_sigma_f_func = nu_sigma_f_ren;
    transitorio.sigma_s_func = params.sigma_s_func;
    transitorio.keff0 = keff0;

    transitorio = transitorio.ensamblar_matrices_cte();

    y0 = transitorio.cond_ini_PC0(Phi0, 0);

    %% Paso 4: Integracion temporal
    paso_t = params.time_delta;
    t = (0 : paso_t : params.t_total)';
    Nt = length(t) - 1;
    K = transitorio.n_precursor;

    y = zeros(length(t), length(y0));
    y(1, :) = y0';
    n_dof = size(transitorio.masas_mat, 1);

    for n = 1 : Nt
        Phi_n = y(n, 1 : G * n_dof)';
        if K > 0
            PC_n = reshape(y(n, G * n_dof + 1 : end), n_dof, K);
        else
            PC_n = [];
        end

        [Phi_act, PC_act] = transitorio.esqPC(paso_t, Phi_n, PC_n, t(n+1));

        if K > 0
            y_act = [Phi_act; PC_act(:)];
        else
            y_act = Phi_act;
        end
        y(n+1, :) = y_act';
    end

    %% Paso 5: Resultados
    transitorio.Phi_t = y(:, 1 : G * n_dof);
    transitorio.Pot_num = transitorio.pot_norm(t, y);

    if bool_graf == true
        transitorio.graficar_t(t, y);
    end
end