%% error_fem_gTrans.m
% Comprueba la consistencia del metodo numerico mediante los bancos de datos
% y grafica la potencia numerica frente a la analitica (si corresponde)

clear; close all; clc;

tol = 1e-3;
metodo = 'fem_2gTransPC';
bool_graf = false;
banco_datos_g = {
    'Mat1D1gTransCte';
    'Mat1D1gTransHomo_absesc';
    'Mat1D1gTransHomo_abslin';
    'Mat1D1gTransHomo_fisesc';
    'Mat1D1gTransHomo_fislin';
    'Mat1D1gTransComp';
    'Mat1D2gTransCte';
    'Mat1D2gTransHomo_absesc';
    'Mat1D2gTransHomo_abslin';
    'Mat1D2gTransHomo_fisesc';
    'Mat1D2gTransHomo_fislin';
    'Mat1D2gTransComp';
    'Mat1D7gTransCte';
};
banco_datos_2g = {
    'Mat1D2gTransCte';
    'Mat1D2gTransHomo_absesc';
    'Mat1D2gTransHomo_abslin';
    'Mat1D2gTransHomo_fisesc';
    'Mat1D2gTransHomo_fislin';
    'Mat1D2gTransComp';
};

check_banco(banco_datos_g, tol, bool_graf, metodo);
fprintf('\n*** FIN ***')