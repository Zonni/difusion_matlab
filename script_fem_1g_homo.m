% DESCRIPCIÓN:
% Script principal para resolver difusión de 1 GRUPO mediante FEM (polinomios orden 3).
% Simula un reactor HOMOGÉNEO de 350 cm (todo el dominio se define como material 2).
% Calcula los 3 primeros k_eff, grafica modos y guarda resultados.
clearvars;
close all; 

L = 350;

% Propiedades físicas (definidas para 2 materiales posibles)
D = [1.446, 0.776];
sigma_a = [0.0077, 0.0244];
nu_sigma_f = [0, 0.0260];

% Configuración de Malla y Materiales
tamano_celdas=[25, 25 + zeros(1, 12), 25];
materiales = [2, 2 + zeros(1, 12), 2]; % Se asigna Material 2 a TODOS los elementos (Homogéneo)
N = 14;
grado_l = 3;

output_file = 'FEM_1G_HOM_N14.mat';

% Inicialización de objetos (Clases de 1 Grupo)
malla = Malla1D(L, N, grado_l, tamano_celdas);
materiales = Materiales1D1g(materiales, D, sigma_a, nu_sigma_f, malla);
elemento = ElementoFinito(grado_l);
problema = ProblemaDifusion1D1g(malla, materiales, elemento);

% Resolución del sistema
problema = problema.ensamblar_matrices_fem();
problema = problema.aplicar_cc();

format long
problema = problema.resolver_autovalor(3);

% Visualización y guardado
problema.graficar([1,2,3]);
save(output_file)
disp(problema.keff)
load('FEM_1G_HOM_N14.mat')
