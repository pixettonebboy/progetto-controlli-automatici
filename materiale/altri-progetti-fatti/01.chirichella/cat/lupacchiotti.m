% PROGETTO CAT 2024/2025 -- Tipologia B traccia 3
% di Filippo Giulietti, Michele Proietti, Renato Eugenio Maria Marziano, Davide Chirichella

clear all; close all; clc;

%% ------------------------------------------------------------------------
%                           FASE DI ANALISI
% -------------------------------------------------------------------------

k = 500;               % elasticità del disco [N/m^2]
beta = 0.5;            % attrito viscoso [N*s/m^2]
alpha = deg2rad(50);   % angolo del giunto di Cardano [rad]
J = 400;               % momento di inerzia della tavola [J*m^2]

% Funzione anonima tau(theta)
tau = @(theta) cos(alpha) ./ (1 - (sin(alpha) .* cos(theta)).^2);


%% DEFINIZIONE FUNZIONI DI STATO
% x = [theta, omega],  u = coppia motore (in ingresso)
% preso dalla traccia, transformato in forma di stato risolto per x1' e x2'
f1 = @(theta, omega, u) omega;
f2 = @(theta, omega, u) (tau(theta) .* u - beta .* omega - k .* theta) / J;


%% EQUILIBRIO DEL SISTEMA
% [rad] angolo di equilibrio
theta_e = deg2rad(100);
x1_e = theta_e;
x2_e = 0; % 0 è la derivata di theta che è di equilibrio
x_e = [ x1_e , x2_e ];
% Calcolo coppia di equilibrio
Cm_e = (k * theta_e) / tau(theta_e);  % corrispondente al u_e


%% LINEARIZZAZIONE. CALCOLO DI A_e E B_e
% Calcolo derivate parziali
delta = 10^(-6); % Incremento per derivate finite

% Matrice A_e: derivate rispetto a [theta, omega]
A_e = zeros(2, 2);
A_e(1, 1) = (f1(theta_e + delta, 0, Cm_e) - f1(theta_e, 0, Cm_e)) / delta; % df1/dtheta
A_e(1, 2) = (f1(theta_e, delta, Cm_e) - f1(theta_e, 0, Cm_e)) / delta;    % df1/domega
A_e(2, 1) = (f2(theta_e + delta, 0, Cm_e) - f2(theta_e, 0, Cm_e)) / delta; % df2/dtheta
A_e(2, 2) = (f2(theta_e, delta, Cm_e) - f2(theta_e, 0, Cm_e)) / delta;    % df2/domega

% Matrice B_e: derivate rispetto a u
B_e = zeros(2, 1);
B_e(1, 1) = (f1(theta_e, 0, Cm_e + delta) - f1(theta_e, 0, Cm_e)) / delta; % df1/du
B_e(2, 1) = (f2(theta_e, 0, Cm_e + delta) - f2(theta_e, 0, Cm_e)) / delta; % df2/du

% Uscita: y = theta
% avendo C = [1, 0] e D = 0
C = [1, 0];
D = 0;

% Costruiamo il modello
sys_ss = ss(A_e, B_e, C, D);




%% ------------------------------------------------------------------------
%                         FASE DI PROGETTAZIONE
% -------------------------------------------------------------------------

%specifiche
e_star = 0.01;     % errore a regime
W = 1.5;           % Ampiezza gradino in riferimento ad errore a regime
D_amp = 1;         % Ampiezza disturbo in riferimento ad errore a regime
Mf_specificato = 33;

T_star = 0.003;    % Tempo di assestamento max
sovr_max = 0.16;   % Sovraelongazione massima (16%)

A_d = 50;          % attenuazione disturbo sull'uscita su banda [ 0, 0.8 ]
Banda_A_d = { 0.001  ;  0.8 };

A_n = 72;          % attenuazione disturbo di misura su banda [1.2 * 10^5, 5 * 10^6]
Banda_A_n = { 1.2 * 10^5  ;  5 * 10^6 };


%% CALCOLO SPECIFICHE DERIVATE
% Coefficiente di smorzamento csi* (dalla sovraelongazione max)
csi_star = abs(log(sovr_max)) / sqrt(pi^2 + (log(sovr_max))^2);
% Margine di fase (in gradi)
Mf = max(csi_star*100, Mf_specificato);


%% DEFINIZIONE ZONE DEL DIAGRAMMA
% Zona proibita legata a Wp
%voglio una certa pulsazione
omegaCmin = 300 / (Mf * T_star);
omega_d_MAX = 0.8;

% Zona proibita legata a d(t)
ZP_d = { [Banda_A_d{1}; Banda_A_d{2}; Banda_A_d{2}; Banda_A_d{1}]  ;  [A_d; A_d; -150; -150] };
% Zona proibita legata a n(t)
ZP_n = { [Banda_A_n{1}; Banda_A_n{2}; Banda_A_n{2}; Banda_A_n{1}] ;  [-A_n; -A_n; 100; 100] };
% Zona proibita legata a omega_c per rispettare la sovraelongazione
ZP_OmegaC = { [Banda_A_d{2}; omegaCmin; omegaCmin; Banda_A_d{2} ]  ;  [0; 0; -150; -150] };
% Zona legata al Margine di Fase
phi_low = -450;      % minimo/massimo sfasamento nel grafico (per funzione di plot)
ZP_Mf = { [omegaCmin; Banda_A_n{1}; Banda_A_n{1}; omegaCmin]  ;  [Mf-180; Mf-180; phi_low; phi_low] };


% Diagramma di Bode di G(s)
s = tf('s');  % variabile di LaPlace
G = tf(sys_ss);  % funzione di trasferimento G del sistema
%disegnaBode(G, 'Bode di G(s)', ZP_d, ZP_n, ZP_OmegaC, ZP_Mf);

%% REGOLATORE STATICO
% 1. Il guadagno statico >= al minimo calcolato per rispettare errore a regime
mu_s_error = (W + D_amp) / e_star;   
mu_s_dist   = 10^(A_d/20);

% Calcolo valori in 0 e Omega MAX
%guadagni dell'impanto
G_0         = abs(evalfr(G, 0));
G_om_d_MAX  = abs(evalfr(G, 1j*omega_d_MAX));

% Definizione regolatore statico
R_s        = max(mu_s_error/G_0, mu_s_dist/G_om_d_MAX); 

% Diagramma esteso
G_e = series(R_s, G);
%disegnaBode(G_e, 'Bode di G_e(s)', ZP_d, ZP_n, ZP_OmegaC, ZP_Mf);

%% REGOLATORE DINAMICO
%un polo e una rete anticipatrice

% 1. Un polo per rispettare attenuazione errori di misura
T1_f = 0.3 * 10^(-4);
R_d_polo = (1)/(1 + T1_f * s);  %NB necessario per attenuare gli errori di misura!

% 2. Rete Anticipatrice (Scenario B visto a lezione)
% - Aumenta la fase per rispettarne il margine
% - Sposta la pulsazione di taglio a frequenze più alte per rispettare omegaCmin

% Calcolo alpha e tau della rete anticipatrice
% Variabile omega_c_star (pulsazione di massimo aumento di fase della rete anticipatrice)
omega_c_star = 1.1 * omegaCmin; % Nuova frequenza centrale

mag_omega_c_star = abs(evalfr(G_e, j * omega_c_star));
arg_omega_c_star = rad2deg(angle(evalfr(G_e, j * omega_c_star)));

M_star = 1 / mag_omega_c_star;
phi_star = Mf + 8 - 180 - arg_omega_c_star;

tau = (M_star - cos(phi_star * pi / 180)) / (omega_c_star * sin(phi_star * pi / 180));
alpha_tau = (cos(phi_star * pi / 180) - 1 / M_star) / (omega_c_star * sin(phi_star * pi / 180));

R_d_antic = (1 + tau * s) / (1 + alpha_tau * s);
%bode(R_d_antic);

% Serie delle parti
R_d = R_d_antic * R_d_polo;

% REGOLATORE COMPLETO
R = series(R_s, R_d);
%disegnaBode(R, 'Bode di R(s)', ZP_d, ZP_n, ZP_OmegaC, ZP_Mf);

%% FUNZIONE DI ANELLO APERTO L
L = series(G, R);
disegnaBode(L, 'Bode di L(s)', ZP_d, ZP_n, ZP_OmegaC, ZP_Mf);

%% ------------------------------------------------------------------------
%            SIMULAZIONI DEL SISTENA IN RISPOSTA AD INGRESSI
% -------------------------------------------------------------------------

F = L/(1 + L); % Funzione di sentitività
S = 1/(1 + L); % Funzione di sentitività complementare

t_sim = 0:0.0000001:0.09;
T_sim = t_sim(end);


%% RISPOSTA AL GRADINO W
[y_w, t_w] = step(W * F, t_sim);

% Plot grafico
figure('Name','Risposta al gradino di riferimento (F)');
plot(t_w, y_w); grid on;
xlabel('Tempo (s)'); ylabel('Ampiezza');
title('Risposta al gradino di riferimento');


%% RISPOSTA AL DISTURBO d(t)
% Costruzione della funzione d(t)
d_t = 0.2 * sin(0.1*1*t_sim) + 0.2 * sin(0.1*2*t_sim) + 0.2 * sin(0.1*3*t_sim) + 0.2 * sin(0.1*4*t_sim);
y_d = lsim(S, d_t, t_sim);

% Plot grafico
figure('Name','Risposta al disturbo d(t)');
plot(t_sim, d_t, 'r--', t_sim, y_d, 'b'); grid on;
legend('d(t)', 'y_d(t)');
xlabel('Tempo (s)'); ylabel('Ampiezza');


%% RISPOSTA AGLI ERRORI DI MISURA n(t)
n_t = 0.2*sin(0.1*1*10^5*t_sim) + 0.2*sin(0.1*2*10^5*t_sim) + 0.2*sin(0.1*3*10^5*t_sim) + 0.2*sin(0.1*4*10^5*t_sim);
y_n = lsim(-F, n_t, t_sim);

% Plot grafico
figure('Name','Risposta agli errori di misura n(t)');
plot(t_sim, n_t, 'r--', t_sim, y_n, 'b'); grid on;
legend('n(t)', 'y_n(t)');
xlabel('Tempo (s)'); ylabel('Ampiezza');


%% RISPOSTA TOTALE con SOVRAPPOSIZIONE
y = y_w + y_d + y_n;
LV = evalfr(W * F, 0);

figure('Name','Risposta totale');
grid on; hold on;
plot(t_sim, y_w, 'r--', t_sim, y, 'b'); grid on;
legend('y_w', 'y(t)');
xlabel('Tempo (s)'); ylabel('Ampiezza');
ylim([0, LV*2]);

% Zona vincolo sovraelongazione
patch([0,T_sim,T_sim,0],[LV*(1+sovr_max), LV*(1+sovr_max), LV*2, LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% Zona vincolo tempo di assestamento al 5%
patch([T_star,T_sim,T_sim,T_star],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_sim,T_sim,T_star],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);
