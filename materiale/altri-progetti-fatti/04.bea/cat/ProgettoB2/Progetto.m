%% MAIN
%% Specifiche globali
close all; clc

% Pulsazione minima e massima nei grafici
omega_plot_min = 1e-4;
omega_plot_max = 1e8;

%% Definizione delle costanti
alpha = deg2rad(25);  % Angolo tra i due alberi
beta = 0.77;          % Coefficiente di attrito viscoso
k = 1000;             % Coefficiente di elasticità del disco
J = 549;              % Momento di inerzia della tavola


%% Specifiche di progetto
% Errore a regime
WW          = 2;
DD          = 1.2;
e_star      = 0.05;

% Attenuazione disturbo sull'uscita
A_d         = 35;
omega_d_min = omega_plot_min;
omega_d_MAX = 0.7;

% Attenuazione disturbo di misura
A_n         = 69;
omega_n_min = 2e5;
omega_n_MAX = 5e6;

% Sovraelongazione massima e tempo d'assestamento al 5%
S_star      = 20;
T_star      = 6e-3;     

% Margine di fase
Mf_esp      = 35;


%% Punto 1
%% Linearizzazione del sistema
Punto1;

%% Punto 2
%% Calcolo e plot della funzione di trasferimento
Punto2;

%% Punto 3
%% Sintesi regolatore da specifiche di progetto
Punto3;

%% Punto 4
%% Test sul sistema sistema linearizzato
Punto4;

