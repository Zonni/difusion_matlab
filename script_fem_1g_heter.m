% DESCRIPCIÓN:
% Script principal para resolver difusión de 1 GRUPO en configuración HETEROGÉNEA.
% Utiliza FEM con polinomios de grado 4.
% Configuración "Reflector-Núcleo-Reflector":
% - Extremos (Celdas 1 y 14): Material 1 (Reflector).
% - Centro (Celdas 2-13): Material 2 (Combustible).
clearvars;
close all; 
%%
L = 350;

% Propiedades Físicas (Mat 1: Reflector, Mat 2: Combustible)
D = [1.446, 0.776];
sigma_a = [0.0077, 0.0244];
nu_sigma_f = [0, 0.0260]; % Mat 1 no tiene fisión

% Definición de la Malla y distribución de materiales
tamano_celdas=[25, 25 + zeros(1, 12), 25];
materiales = [1, 2 + zeros(1, 12), 1]; % Mat 1 en extremos, Mat 2 en el centro
N = 14;
grado_l = 4;

output_file = 'FEM_1G_HET_N14.mat';

% Inicialización de objetos
malla = Malla1D(L, N, grado_l, tamano_celdas);
materiales = Materiales1D1g(materiales, D, sigma_a, nu_sigma_f, malla);
elemento = ElementoFinito(grado_l);
problema = ProblemaDifusion1D1g(malla, materiales, elemento);

% Resolución del sistema
problema = problema.ensamblar_matrices_fem();
problema = problema.aplicar_cc();

format long
problema = problema.resolver_autovalor(3); % Calcular 3 primeros modos

% Resultados
problema.graficar([1,2,3]);
save(output_file)
disp(problema.keff)
load('FEM_1G_HET_N14.mat')
