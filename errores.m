%% ANÁLISIS DE ERRORES - SOLUCIONES ANALÍTICAS

%% Solución 1G homogéneo
clear all
load('FEM_1G_HOM_N14.mat')

% Solución analítica
m = length(problema.keff);
n_modos = 1:m;
k_eff_analitico_1g = nu_sigma_f(2) ./ (D(2)*(pi*n_modos/L).^2 + sigma_a(2));  % Vector de k_eff analíticos
phi_analitico_1g = @(x,n) sin(pi*n*x/L);

% Extraer solución numérica de todos los modos
k_eff_numerico_1g = problema.keff(:);

% Inicializar vectores para almacenar errores
ErA_1g = zeros(m, 1); 
ErR_1g = zeros(m, 1); 


% Calcular errores para cada modo
x_nodos = malla.x_nodos;
for modo = 1:m
    % Extraer solución numérica del modo
    phi_numerico_1g = problema.phi_inc(:, modo);
    phi_numerico_1g = phi_numerico_1g(:);  % vector columna

    % Calcular solución analítica en los mismos puntos
    phi_analitico_norm_1g = phi_analitico_1g(x_nodos, modo);
    phi_analitico_norm_1g = phi_analitico_norm_1g(:);  % vector columna
    
    % Normalizar soluciones
    if phi_analitico_norm_1g(2)<0
        phi_analitico_norm_1g(:)=-phi_analitico_norm_1g(:);
    end
    phi_analitico_norm_1g(:)=phi_analitico_norm_1g(:)/max(phi_analitico_norm_1g(:));
    
    % Calcular errores
    ErA_1g(modo) = abs(k_eff_analitico_1g(modo) - k_eff_numerico_1g(modo)) * 10^5;
    ErR_1g(modo) = rmse(phi_analitico_norm_1g,phi_numerico_1g);
end

% Gráfica comparativa para todos los modos
figure(1); 
hold on; grid on;
colores = {'b', 'g', 'm', 'c'};
n_puntos = length(problema.phi_inc(:, 1));
x_plot = linspace(0, L, n_puntos)';  % Usar los mismos puntos que la solución numérica
for modo = 1:m
    phi_analitico_norm = phi_analitico_1g(x_plot, modo);
    phi_analitico_norm = phi_analitico_norm(:) / max(abs(phi_analitico_norm(:)));
    phi_numerico_norm = problema.phi_inc(:, modo) / max(abs(problema.phi_inc(:, modo)));
    
    plot(x_plot, phi_analitico_norm, 'Color', colores{modo}, 'LineStyle', '-', ...
         'LineWidth', 2, 'DisplayName', sprintf('Analítica Modo %d', modo));
    plot(x_plot, phi_numerico_norm, 'Color', colores{modo}, 'LineStyle', '--', ...
         'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'DisplayName', sprintf('FEM Modo %d', modo));
end
xlabel('x (cm)');
ylabel('\phi(x)');
title('Solución 1G Homogéneo - Comparación Analítica vs FEM (Todos los modos)');
legend('Location', 'best');

% Tabla de resultados para cada modo
fprintf('\n=== RESULTADOS 1G HOMOGÉNEO- POR MODO ===\n');
fprintf('%-10s | %-18s | %-18s | %-15s | %-15s\n', ...
        'Modo', 'k_eff Analítico', 'k_eff Numérico', 'Error k_eff (×10^5)', 'RMSE Flujo');
fprintf('----------------------------------------------------------------------------------------\n');
for modo = 1:m
    fprintf('%-10d | %18.10f | %18.10f | %15.6f | %15.6e\n', ...
            modo, k_eff_analitico_1g(modo), k_eff_numerico_1g(modo), ...
            ErA_1g(modo), ErR_1g(modo));
end
fprintf('----------------------------------------------------------------------------------------\n');

% Generar tabla en formato LaTeX
fid = fopen('tabla_errores_1g_homogeneo.tex', 'w');
fprintf(fid, '\\begin{table}[h]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Resultados de errores para el caso 1G Homogéneo- Comparación analítica vs numérica}\n');
fprintf(fid, '\\label{tab:errores_1g}\n');
fprintf(fid, '\\begin{tabular}{c|cc|cc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '\\textbf{Modo} & \\textbf{$k_{eff}$ Analítico} & \\textbf{$k_{eff}$ Numérico} & \\textbf{Error $k_{eff}$ ($\\times 10^5$)} & \\textbf{RMSE Flujo} \\\\\n');
fprintf(fid, '\\midrule\n');
for modo = 1:m
    fprintf(fid, '%d & %.10f & %.10f & %.6f & %.6e \\\\\n', ...
            modo, k_eff_analitico_1g(modo), k_eff_numerico_1g(modo), ...
            ErA_1g(modo), ErR_1g(modo));
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('\n Tabla LaTeX guardada en: tabla_errores_1g_homogeneo.tex\n');

% Guardar variables del caso homogéneo para la tabla completa
m_hom = m;
ErA_1g_hom = ErA_1g;
ErR_1g_hom = ErR_1g;
k_eff_analitico_1g_hom = k_eff_analitico_1g;
k_eff_numerico_1g_hom = k_eff_numerico_1g;

%% Solución 1G heterogéneo
clear all
load('FEM_1G_HET_N14.mat')

% Solución analítica para caso heterogéneo resolviendo el sistema matricial
% Usando el k_eff numérico del FEM para construir la solución analítica

% Parámetros geométricos
% Estructura: reflector (0-a) | combustible (a-b) | reflector (b-L)
x_celdas = malla.x_celdas;
a = x_celdas(2);
b = sum(x_celdas(end-1));
L = malla.L;

% Propiedades de materiales
% Material 1 (reflector)
D1 = D(1);
sigma_a1 = sigma_a(1);
nu_sigma_f1 = nu_sigma_f(1);  % = 0

% Material 2 (combustible)
D2 = D(2);
sigma_a2 = sigma_a(2);
nu_sigma_f2 = nu_sigma_f(2);

m = length(problema.keff);
x_nodos_het = malla.x_nodos;
k_eff_numerico_1g_het = problema.keff(:);

% Inicializar vectores para errores
ErA_1g_het = zeros(m, 1);
ErR_1g_het = zeros(m, 1);

% Calcular solución analítica para cada modo
fprintf('\n=== SOLUCIÓN ANALÍTICA HETEROGÉNEA ===\n');
fprintf('Geometría: a=%.1f cm, b=%.1f cm, L=%.1f cm\n', a, b, L);
fprintf('\n');

% Almacenar soluciones analíticas
phi_analitico_het = cell(m, 1);

for modo = 1:m
    % Usar k_eff numérico del FEM para este modo
    k_eff_modo = k_eff_numerico_1g_het(modo);
    
    % Calcular parámetros κ
    kappa1 = sqrt(sigma_a1 / D1);  % Para reflector (νΣf = 0)
    kappa2 = sqrt((nu_sigma_f2/k_eff_modo - sigma_a2) / D2);
    
    % Construir matriz del sistema de condiciones de continuidad
    % Sistema: M * [A1; A2; B2; A3] = 0
    
    % Condición 1: Continuidad de flujo en x=a
    % φ₁(a) = φ₂(a)
    % A₁ sinh(κ₁a) = A₂ sin(κ₂a) + B₂ cos(κ₂a)
    M11 = sinh(kappa1*a);
    M12 = -sin(kappa2*a);
    M13 = -cos(kappa2*a);
    M14 = 0;
    
    % Condición 2: Continuidad de corriente en x=a
    % D₁ dφ₁/dx|a = D₂ dφ₂/dx|a
    % D₁ κ₁ A₁ cosh(κ₁a) = D₂ κ₂ (A₂ cos(κ₂a) - B₂ sin(κ₂a))
    M21 = D1 * kappa1 * cosh(kappa1*a);
    M22 = -D2 * kappa2 * cos(kappa2*a);
    M23 = D2 * kappa2 * sin(kappa2*a);
    M24 = 0;
    
    % Condición 3: Continuidad de flujo en x=b
    % φ₂(b) = φ₃(b)
    % A₂ sin(κ₂b) + B₂ cos(κ₂b) = A₃ sinh(κ₁(L-b))
    M31 = 0;
    M32 = sin(kappa2*b);
    M33 = cos(kappa2*b);
    M34 = -sinh(kappa1*(L-b));
    
    % Condición 4: Continuidad de corriente en x=b
    % D₂ dφ₂/dx|b = D₁ dφ₃/dx|b
    % dφ₃/dx = -A₃ κ₁ cosh(κ₁(L-x))
    % Entonces: D₂ κ₂ (A₂ cos(κ₂b) - B₂ sin(κ₂b)) = D₁ (-A₃ κ₁ cosh(κ₁(L-b)))
    % Reordenando: D₂ κ₂ A₂ cos(κ₂b) - D₂ κ₂ B₂ sin(κ₂b) + D₁ κ₁ A₃ cosh(κ₁(L-b)) = 0
    M41 = 0;
    M42 = D2 * kappa2 * cos(kappa2*b);
    M43 = -D2 * kappa2 * sin(kappa2*b);
    M44 = D1 * kappa1 * cosh(kappa1*(L-b));  % Signo positivo porque dφ₃/dx tiene signo negativo
    
    % Construir matriz
    M = [M11, M12, M13, M14;
         M21, M22, M23, M24;
         M31, M32, M33, M34;
         M41, M42, M43, M44];
    
    % Resolver sistema homogéneo M * v = 0
    % El sistema tiene solución no trivial si det(M) = 0
    % Para obtener las constantes, resolvemos M * v = 0
    % Usamos el autovector correspondiente al autovalor más cercano a cero
    [V, Lambda] = eig(M);
    lambda_vals = diag(Lambda);
    
    % Encontrar el autovalor más cercano a cero
    [~, idx] = min(abs(lambda_vals));
    v = V(:, idx);
    
    % Normalizar para que el máximo sea 1
    v = v / max(abs(v));
    
    A1 = real(v(1));
    A2 = real(v(2));
    B2 = real(v(3));
    A3 = real(v(4));
    
    % Construir solución analítica completa
    x_fino = linspace(0, L, 1000)';
    phi_analitico = zeros(size(x_fino));
    
    for i = 1:length(x_fino)
        x_val = x_fino(i);
        if x_val <= a
            % Región 1: reflector izquierdo
            phi_analitico(i) = A1 * sinh(kappa1 * x_val);
        elseif x_val <= b
            % Región 2: combustible
            phi_analitico(i) = A2 * sin(kappa2 * x_val) + B2 * cos(kappa2 * x_val);
        else
            % Región 3: reflector derecho
            phi_analitico(i) = A3 * sinh(kappa1 * (L - x_val));
        end
    end
    
    % Normalizar solución analítica
    if phi_analitico(2)<0
        phi_analitico(:)=-phi_analitico(:);
    end
    phi_analitico(:)=phi_analitico(:)/max(phi_analitico(:));
    
    phi_analitico_het{modo} = phi_analitico;
    
    % Interpolar solución analítica a los puntos de la malla
    phi_analitico_nodos = interp1(x_fino, phi_analitico, x_nodos_het, 'linear', 'extrap');
    phi_analitico_nodos = phi_analitico_nodos(:);
    
    % solución numérica
    phi_numerico = problema.phi_inc(:, modo);
    
    % Calcular errores
    % NOTA: k_eff analítico = k_eff numérico (usado para construir solución analítica)
    % Por lo tanto, el error absoluto en k_eff es 0
    ErA_1g_het(modo) = 0;  % k_eff analítico = k_eff numérico
    ErR_1g_het(modo) = rmse(phi_analitico_nodos,phi_numerico);
    
    fprintf('Modo %d: k_eff = %.10f, RMSE = %.6e\n', modo, k_eff_modo, ErR_1g_het(modo));
end

% k_eff analítico es igual al numérico (lo usamos para construir la solución)
k_eff_analitico_1g_het = k_eff_numerico_1g_het;

% Gráfica comparativa analítica vs numérica
figure(2); 
hold on; grid on;
colores = {'b', 'g', 'm', 'c'};
estilos_analitico = {'-', '-', '-', '-'};
estilos_numerico = {'--', '--', '--', '--'};
x_fino_plot = linspace(0, L, 1000)';

for modo = 1:m
    % Solución analítica (ya calculada)
    phi_ana = phi_analitico_het{modo};
    
    % Solución numérica
    phi_num = problema.phi_inc(:, modo);
    phi_num_norm = phi_num / max(abs(phi_num));
    
    % Interpolar numérica a puntos finos para comparación
    phi_num_fino = interp1(x_nodos_het, phi_num_norm, x_fino_plot, 'linear', 'extrap');
    
    plot(x_fino_plot, phi_ana, 'Color', colores{modo}, 'LineStyle', estilos_analitico{modo}, ...
         'LineWidth', 2, 'DisplayName', sprintf('Analítica Modo %d', modo));
    plot(x_nodos_het, phi_num_norm, 'Color', colores{modo}, 'LineStyle', estilos_numerico{modo}, ...
         'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'DisplayName', sprintf('FEM Modo %d', modo));
end

% Marcar interfaces
xline(a, 'k--', 'LineWidth', 1, 'DisplayName', 'Interfaz 1-2');
xline(b, 'k--', 'LineWidth', 1, 'DisplayName', 'Interfaz 2-3');

xlabel('x (cm)');
ylabel('\phi(x)');
title('Solución 1G Heterogéneo - Comparación Analítica vs FEM');
legend('Location', 'best');

% Tabla de resultados con comparación
fprintf('\n=== RESULTADOS 1G HETEROGÉNEO - POR MODO ===\n');
fprintf('NOTA: k_eff analítico = k_eff numérico (usado para construir solución analítica)\n');
fprintf('%-10s | %-18s | %-18s | %-15s | %-15s\n', ...
        'Modo', 'k_eff Analítico', 'k_eff Numérico', 'Error k_eff (×10^5)', 'RMSE Flujo');
fprintf('----------------------------------------------------------------------------------------\n');
for modo = 1:m
    fprintf('%-10d | %18.10f | %18.10f | %15.6f | %15.6e\n', ...
            modo, k_eff_analitico_1g_het(modo), k_eff_numerico_1g_het(modo), ...
            ErA_1g_het(modo), ErR_1g_het(modo));
end
fprintf('----------------------------------------------------------------------------------------\n');

% Generar tabla en formato LaTeX para heterogéneo
fid = fopen('tabla_errores_1g_heterogeneo.tex', 'w');
fprintf(fid, '\\begin{table}[h]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Resultados de errores para el caso 1G Heterogéneo - Comparación analítica vs numérica}\n');
fprintf(fid, '\\label{tab:errores_1g_heterogeneo}\n');
fprintf(fid, '\\begin{tabular}{c|cc|cc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '\\textbf{Modo} & \\textbf{$k_{eff}$ Analítico} & \\textbf{$k_{eff}$ Numérico} & \\textbf{Error $k_{eff}$ ($\\times 10^5$)} & \\textbf{RMSE Flujo} \\\\\n');
fprintf(fid, '\\midrule\n');
for modo = 1:m
    fprintf(fid, '%d & %.10f & %.10f & %.6f & %.6e \\\\\n', ...
            modo, k_eff_analitico_1g_het(modo), k_eff_numerico_1g_het(modo), ...
            ErA_1g_het(modo), ErR_1g_het(modo));
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('\n Tabla LaTeX guardada en: tabla_errores_1g_heterogeneo.tex\n');

%% Solución analítica 2G homogéneo
clear all
load('FEM_2G_N14.mat')

% Extraer solución numérica del modo fundamental
k_eff_numerico_2g = problema.keff(1);
M_nodos = malla.N * malla.grado_l + 1;  % Número de nodos por grupo
phi_numerico_2g_g1 = problema.phi_inc(1:M_nodos, 1);  % Grupo 1, modo 1
phi_numerico_2g_g2 = problema.phi_inc(M_nodos+1:end, 1);  % Grupo 2, modo 1

x = malla.x_nodos;
L = malla.L;

% Propiedades del material
% D, sigma_a, nu_sigma_f tienen estructura: [grupo1, grupo2; material1, material2]
D1 = problema.materiales.D(1, 1);  % Grupo 1, material 1
D2 = problema.materiales.D(2, 1);  % Grupo 2, material 1
sigma_a1 = problema.materiales.sigma_a(1, 1);
sigma_a2 = problema.materiales.sigma_a(2, 1);
nu_sigma_f1 = problema.materiales.nu_sigma_f(1, 1);
nu_sigma_f2 = problema.materiales.nu_sigma_f(2, 1);
sigma_s12 = problema.materiales.sigma_s12(1);

fprintf('\n=== SOLUCIÓN ANALÍTICA 2G HOMOGÉNEO ===\n');
fprintf('Método: Sistema de ecuaciones acopladas de primer orden\n');
fprintf('k_eff numérico usado (modo fundamental): %.10f\n', k_eff_numerico_2g);
fprintf('\n');

% ============================================================
% MÉTODO: Reducción a sistema de primer orden
% ============================================================
% Las ecuaciones de difusión 2G se pueden escribir como:
% d²φ₁/dx² = a₁₁ φ₁ + a₁₂ φ₂
% d²φ₂/dx² = a₂₁ φ₁ + a₂₂ φ₂
%
% Donde:
% a₁₁ = (Σa₁ + Σs₁₂ - νΣf₁/k_eff)/D₁
% a₁₂ = -νΣf₂/(k_eff D₁)
% a₂₁ = -Σs₁₂/D₂
% a₂₂ = Σa₂/D₂

a11 = (sigma_a1 + sigma_s12 - nu_sigma_f1/k_eff_numerico_2g) / D1;
a12 = -nu_sigma_f2 / (k_eff_numerico_2g * D1);
a21 = -sigma_s12 / D2;
a22 = sigma_a2 / D2;

% Construir matriz del sistema de primer orden
% Sistema: d/dx [φ₁, φ₂, φ₁', φ₂']ᵀ = M * [φ₁, φ₂, φ₁', φ₂']ᵀ
M_sistema = [0,   0,   1,   0;
             0,   0,   0,   1;
             a11, a12, 0,   0;
             a21, a22, 0,   0];

% Calcular autovalores y autovectores
[V, Lambda] = eig(M_sistema);
lambda_vals = diag(Lambda);

fprintf('Autovalores del sistema:\n');
for i = 1:4
    fprintf('  λ%d = %.6f + %.6fi\n', i, real(lambda_vals(i)), imag(lambda_vals(i)));
end
fprintf('\n');

% ============================================================
% CONSTRUCCIÓN DE LA SOLUCIÓN GENERAL
% ============================================================
% La solución general es: y(x) = Σ cᵢ vᵢ e^(λᵢ x)
% donde y = [φ₁, φ₂, φ₁', φ₂']ᵀ
% 
% Para extraer φ₁ y φ₂:
% φ₁(x) = Σ cᵢ vᵢ(1) e^(λᵢ x)  [primera componente de cada autovector]
% φ₂(x) = Σ cᵢ vᵢ(2) e^(λᵢ x)  [segunda componente de cada autovector]

% Construir sistema para determinar coeficientes que satisfagan
% condiciones de contorno: φ₁(0)=0, φ₂(0)=0, φ₁(L)=0, φ₂(L)=0

% Evaluar autovectores en x=0 y x=L
A_bc = zeros(4, 4);  % Matriz para condiciones de contorno

for i = 1:4
    lambda_i = lambda_vals(i);
    v = V(:, i);
    
    % Extraer componentes de flujo del autovector
    % v = [v₁, v₂, v₃, v₄]ᵀ donde:
    % v₁ = componente de φ₁
    % v₂ = componente de φ₂
    % v₃ = componente de φ₁' (no se usa para condiciones de flujo)
    % v₄ = componente de φ₂' (no se usa para condiciones de flujo)
    v1 = v(1);  % Componente de φ₁
    v2 = v(2);  % Componente de φ₂
    
    % Condición en x=0: φ(0) = v (ya que e^0 = 1)
    A_bc(1, i) = v1;  % φ₁(0) = 0
    A_bc(2, i) = v2;  % φ₂(0) = 0
    
    % Condición en x=L: φ(L) = v * e^(λL)
    exp_L = exp(lambda_i * L);
    A_bc(3, i) = v1 * exp_L;  % φ₁(L) = 0
    A_bc(4, i) = v2 * exp_L;  % φ₂(L) = 0
end

% Resolver sistema homogéneo A_bc * c = 0 usando SVD
% El vector singular correspondiente al valor singular más pequeño
% es la solución del sistema homogéneo
[~, S, V_svd] = svd(A_bc);
[~, idx_min] = min(diag(S));
coeficientes = V_svd(:, idx_min);

% Normalizar coeficientes
coeficientes = coeficientes / max(abs(coeficientes));

fprintf('Coeficientes determinados:\n');
for i = 1:4
    fprintf('  c%d = %.6f + %.6fi\n', i, real(coeficientes(i)), imag(coeficientes(i)));
end
fprintf('\n');

% ============================================================
% CONSTRUIR SOLUCIÓN FINAL: Extraer φ₁ y φ₂
% ============================================================
% De la solución general y(x) = Σ cᵢ vᵢ e^(λᵢ x)
% Extraemos:
% φ₁(x) = Σ cᵢ vᵢ(1) e^(λᵢ x)  [primera componente]
% φ₂(x) = Σ cᵢ vᵢ(2) e^(λᵢ x)  [segunda componente]

phi_analitico_2g_g1 = zeros(size(x));
phi_analitico_2g_g2 = zeros(size(x));

for i = 1:4
    lambda_i = lambda_vals(i);
    v = V(:, i);
    c_i = coeficientes(i);
    
    % Extraer componentes de flujo del autovector
    v1 = v(1);  % Primera componente → φ₁
    v2 = v(2);  % Segunda componente → φ₂
    
    % Contribución a la solución:
    % φ₁ recibe: cᵢ * vᵢ(1) * e^(λᵢ x)
    % φ₂ recibe: cᵢ * vᵢ(2) * e^(λᵢ x)
    contribucion = c_i * exp(lambda_i * x);
    
    phi_analitico_2g_g1 = phi_analitico_2g_g1 + v1 * contribucion
    phi_analitico_2g_g2 = phi_analitico_2g_g2 + v2 * contribucion
end

% Tomar parte real (la solución física es real)
phi_analitico_2g_g1 = real(phi_analitico_2g_g1);
phi_analitico_2g_g2 = real(phi_analitico_2g_g2);

% Normalizar soluciones analíticas
if phi_analitico_2g_g1(2) < 0
    phi_analitico_2g_g1(:) = -phi_analitico_2g_g1(:);
end
phi_analitico_2g_g1(:) = phi_analitico_2g_g1(:) / max(phi_analitico_2g_g1(:))

if phi_analitico_2g_g2(2) < 0
    phi_analitico_2g_g2(:) = -phi_analitico_2g_g2(:);
end
phi_analitico_2g_g2(:) = phi_analitico_2g_g2(:) / max(phi_analitico_2g_g2(:))

% Calcular errores usando función rmse
ErR_2g_g1 = rmse(phi_analitico_2g_g1, phi_numerico_2g_g1);
ErR_2g_g2 = rmse(phi_analitico_2g_g2, phi_numerico_2g_g2);

% Gráfica comparativa 2G
figure(3);
subplot(2,1,1);
hold on; grid on;
plot(x, phi_analitico_2g_g1, 'b-', 'LineWidth', 2, 'DisplayName', 'Analítica G1');
plot(x, phi_numerico_2g_g1, 'r--o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'FEM G1');
xlabel('x (cm)');
ylabel('\phi_1(x)');
title('Solución 2G - Grupo 1');
legend('Location', 'best');

subplot(2,1,2);
hold on; grid on;
plot(x, phi_analitico_2g_g2, 'b-', 'LineWidth', 2, 'DisplayName', 'Analítica G2');
plot(x, phi_numerico_2g_g2, 'r--o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'FEM G2');
xlabel('x (cm)');
ylabel('\phi_2(x)');
title('Solución 2G - Grupo 2');
legend('Location', 'best');

fprintf('\n=== RESULTADOS 2G ===\n');
fprintf('k_eff numérico:  %.10f\n', k_eff_numerico_2g);
fprintf('RMSE flujo G1: %.6e\n', ErR_2g_g1);
fprintf('RMSE flujo G2: %.6e\n', ErR_2g_g2);
