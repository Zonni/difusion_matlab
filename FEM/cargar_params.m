function params = cargar_params(banco)
%   Carga los parametros del banco de datos especificado

    % Solicitar opcion mediante el nombre del archivo
    nombre_archivo = banco;
    params = feval(nombre_archivo);

    % Mostrar un resumen de los parametros cargados
    fprintf('\n--- Banco de datos "%s" ---\n', nombre_archivo);
    fprintf('Longitud del reactor: %.1f cm\n', params.L);
    fprintf('Numero de celdas: %d\n', params.N);
    fprintf('Longitud de las celdas: %s\n', mat2str(params.tam_celdas));
    fprintf('Materiales por celda: %s\n', mat2str(params.materiales));
    fprintf('Tiempo total: %.2f s\n', params.t_total);
    if isempty(params.beta)
        fprintf('Sin precursores\n');
    else
        fprintf('%d grupos de precursores.\n', length(params.beta));
    end
end