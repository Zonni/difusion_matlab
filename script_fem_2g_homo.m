% DESCRIPCIÓN:
% Script principal para resolver el problema de difusión de 2 GRUPOS con FEM (orden 3).
% Configuración HOMOGÉNEA:
% - Todo el dominio (elementos 1 a N) se define como Material 1 (Núcleo).

clearvars;
close all; 

L = 350;
N = 14;
grado_l = 3;

% Definición de la geometría y distribución de materiales
tamano_celdas = [25, 25 + zeros(1, 12), 25];

%Todo el reactor es Material 1
materiales = ones(1, N); 

% Propiedades Físicas (dejamos las del caso heterogéno)
% (Filas: Grupo 1, 2 | Columnas: Material 1, 2)
D         = [1.50,    2.0    % G1 
            0.40,    0.3];   % G2 
sigma_a    = [0.010,  0.00  
    0.085,  0.01];
sigma_s12    = [0.020,  0.04]; % Scattering 1->2
nu_sigma_f = [0.000,   0.0   %  G1
 0.135,   0.0];              % G2 

output_file = 'FEM_2G_HOM_N14.mat';

% Inicialización de objetos
malla = Malla1D(L, N, grado_l, tamano_celdas);
materiales = Materiales1D2g(materiales, D, sigma_a, nu_sigma_f, sigma_s12, malla);
elemento = ElementoFinito(grado_l);
problema = ProblemaDifusion1D2g(malla, materiales, elemento);

% Resolución
problema = problema.ensamblar_matrices_fem();
problema = problema.aplicar_cc();
problema = problema.resolver_autovalor(3);

% Resultados
problema.graficar([1,2,3]);
save(output_file)
disp(problema.keff)
load('FEM_2G_HOM_N14.mat')