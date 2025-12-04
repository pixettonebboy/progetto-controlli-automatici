%% =========================================================================
%  PROGETTO CAT 2024/2025 -- Tipologia B, traccia 1
%  Punto 1: Modellazione non lineare, punto di equilibrio, linearizzazione
%           e analisi stabilità sistema linearizzato
% =========================================================================

% punto 1 --> {10 => 150}
% punto 2 --> {150 => 360}
% punto 3 --> {360 => 850}
% punto 4 --> 850 => 1000}

clear; clc; close all;


%% ========================================================================
%  PARTE 1 - PARAMETRI, MODELLO, EQUILIBRIO E LINEARIZZAZIONE
% ========================================================================

%% ------------------------------------------------------------------------
%  PARAMETRI FISICI DEL SISTEMA
% ------------------------------------------------------------------------
% Stato:
%   x = [TR; Tout]
%       TR   : temperatura resistenza (°C)
%       Tout : temperatura aria in uscita (°C)
%
% Ingresso:
%   u = PE   : potenza elettrica (W)
%
% Uscita:
%   y = Tout

hR     = 50;          % coefficiente di convezione (W/(m^2 °C))
AR     = 0.07;        % area di scambio termico (m^2)
cR     = 840.8;       % capacità termica specifica resistenza (J/(kg °C))
cA     = 1010;        % capacità termica specifica aria (J/(kg °C))
mR     = 2.542;       % massa resistenza (kg)
mA     = 0.1041;      % massa aria (kg)
m_dotA = 0.2;         % portata massica aria (kg/s)
Tin    = 28;          % temperatura aria ingresso (°C)
kappa  = 3e-3;        % coefficiente di variazione della resistenza con la temperatura (1/°C)

% Coefficienti utili
aR   = hR * AR;       % coeff. scambio termico globale (W/°C)
den1 = mR * cR;       % capacità termica resistenza (J/°C)
den2 = mA * cA;       % capacità termica aria (J/°C)


%% ------------------------------------------------------------------------
%  MODELLO NON LINEARE (equazioni di stato)
% ------------------------------------------------------------------------
%   den1 * dTR/dt   = aR (Tout - TR) + u/(1 + kappa*TR)
%   den2 * dTout/dt = m_dotA*cA(Tin - Tout) + aR (TR - Tout)

f1 = @(TR, Tout, u) ( aR*(Tout - TR) + u./(1 + kappa*TR) ) / den1;
f2 = @(TR, Tout, u) ( m_dotA*cA*(Tin - Tout) + aR*(TR - Tout) ) / den2;


%% ------------------------------------------------------------------------
%  PUNTO DI EQUILIBRIO (dalla traccia)
% ------------------------------------------------------------------------

% Punto di equilibrio (dato dalla traccia / tabella)
TR_e   = 175;          % [C] temperatura elemento riscaldante all'equilibrio
Tout_e = 80;           % [C] temperatura aria in uscita all'equilibrio
x_e    = [TR_e; Tout_e];

% Ingresso di equilibrio u_e: si impone dTR/dt = 0 all'equilibrio
%   0 = aR*(Tout_e - TR_e) + u_e/(1 + kappa*TR_e)
KTR_e = 1 + kappa * TR_e;
u_e   = - aR * (Tout_e - TR_e) * KTR_e;   % [W] 

% Verifica nel modello non lineare
f1_eq = f1(TR_e, Tout_e, u_e);
f2_eq = f2(TR_e, Tout_e, u_e);

%% ------------------------------------------------------------------------
%  LINEARIZZAZIONE 
% ------------------------------------------------------------------------

% Derivate analitiche di f1
df1_dTR   = (-aR - u_e * kappa / KTR_e^2) / den1;
df1_dTout =  aR / den1;
df1_du    =  1 / (KTR_e * den1);

% Derivate analitiche di f2
df2_dTR   =  aR / den2;
df2_dTout = -(m_dotA*cA + aR) / den2;
df2_du    =  0;

delta = 1e-3;          % passo per l'approssimazione delle derivate

A = zeros(2,2);
B = zeros(2,1);

% df1/dTR
A(1,1) = ( f1(TR_e + delta, Tout_e, u_e) - f1(TR_e - delta, Tout_e, u_e) )/ (2*delta);
% df1/dTout
A(1,2) = ( f1(TR_e, Tout_e + delta, u_e) - f1(TR_e, Tout_e - delta, u_e) ) / (2*delta);

% df2/dTR
A(2,1) = ( f2(TR_e + delta, Tout_e, u_e) - f2(TR_e - delta, Tout_e, u_e) ) / (2*delta);
% df2/dTout
A(2,2) = ( f2(TR_e, Tout_e + delta, u_e) - f2(TR_e, Tout_e - delta, u_e) )/ (2*delta);

% Derivate rispetto all'ingresso u
% df1/du
B(1,1) = ( f1(TR_e, Tout_e, u_e + delta) - f1(TR_e, Tout_e, u_e - delta) )/ (2*delta);
% df2/du
B(2,1) = ( f2(TR_e, Tout_e, u_e + delta) - f2(TR_e, Tout_e, u_e - delta) )/ (2*delta);

% Matrici di uscita (si misura solo Tout)
C = [0 1];
D = 0;

% Sistema linearizzato e poli
lambda = eig(A);


%% ========================================================================
%  PARTE 2 - OUTPUT 
% ========================================================================

fprintf('\n============================================================\n');
fprintf('              PUNTO 1 - RISULTATI PRINCIPALI\n');
fprintf('============================================================\n\n');

fprintf('PUNTO DI EQUILIBRIO:\n');
fprintf('  TR_e   = %8.3f  [°C]\n', TR_e);
fprintf('  Tout_e = %8.3f  [°C]\n', Tout_e);
fprintf('  u_e    = %8.3f  [W]\n\n', u_e);

fprintf('VERIFICA EQUILIBRIO:\n');
fprintf('  f1(TR_e, Tout_e, u_e) = %+10.3e  [°C/s]\n', f1_eq);
fprintf('  f2(TR_e, Tout_e, u_e) = %+10.3e  [°C/s]\n', f2_eq);
if max(abs(f1_eq), abs(f2_eq)) > 1e-6
    fprintf('Le equazioni di stato non soddisfano le condizioni di equilibrio.\n');
end


fprintf('MATRICE A:\n'); disp(A);
fprintf('MATRICE B:\n'); disp(B);
fprintf('MATRICE C:\n'); disp(C);
fprintf('MATRICE D:\n'); disp(D);

fprintf('AUTOVALORI DEL SISTEMA LINEARIZZATO:\n');
disp(lambda);

if all(real(lambda) < 0)
    fprintf('Gli autovalori hanno parte reale negativa: il sistema linearizzato è asintoticamente stabile.\n');
else
    fprintf('Alcuni autovalori presentano parte reale non negativa: il sistema linearizzato non è asintoticamente stabile.\n');
end


%% ======================= PUNTO 2 ==========================
%   FUNZIONE DI TRASFERIMENTO G(s) E ANALISI IN FREQUENZA
% ===========================================================

%% ========================================================================
%  PARTE 1 - CALCOLI E PARAMETRI
% ========================================================================

% Dal modello di stato (A,B,C,D) si ricava G(s) = ΔT_out(s) / ΔPE(s)

% Spazio di stato -> funzione di trasferimento
s=tf('s');
[N,D]=ss2tf(A,B,C,D);
G=tf(N,D);
zpk(G)

% Guadagno statico e poli/zeri
G0 = dcgain(G);      % guadagno statico G(0)
p  = pole(G);        % poli
z  = zero(G);        % zeri

% ----------------- Specifiche del progetto -----------------
Mf_min   = 50;       % [deg] margine di fase minimo richiesto
T_star   = 0.01;     % [s]   tempo di assestamento massimo (5%)
Sovr_max = 0.11;     % [-]   sovraelongazione massima (11%)

A_d      = 50;       % [dB] attenuazione disturbo 
A_n      = 60;       % [dB] attenuazione rumore 

omega_dist_max = 0.4;     % [rad/s] limite superiore banda disturbo
omega_n_min    = 8e4;     % [rad/s] limite inferiore banda rumore
omega_n_max    = 9e6;     % [rad/s] limite superiore banda rumore

% % Smorzamento equivalente ricavato dalla sovraelongazione percentuale
csi_star = abs(log(Sovr_max)) / sqrt(pi^2 + (log(Sovr_max))^2);

% Margine di fase desiderato (in gradi)
Mf_des = max(csi_star * 100, Mf_min);

% Stima euristica pulsazione di taglio minima (da Ta,5 ≈ 3/(ξ ω_n))
omega_c_min = 300 / (Mf_des * T_star);  % [rad/s]

% ----------------- Preparazione dati Bode ------------------

omega_plot_min = 1e-3;
omega_plot_max = 1e7;

[MagG, PhaseG, wG] = bode(G, {omega_plot_min, omega_plot_max});
MagG       = squeeze(MagG);
PhaseG_deg = squeeze(PhaseG);          % [deg]
MagG_dB    = 20*log10(MagG);

% Limiti per gli assi (modulo e fase) basati sui dati
mag_min = min(MagG_dB) - 20;
mag_max = max(MagG_dB) + 20;
phi_min = min(PhaseG_deg) - 20;
phi_max = max(PhaseG_deg) + 20;


%% ========================================================================
%  PARTE 2 - OUTPUT TESTUALE E GRAFICA
% ========================================================================

%% ----------------- REPORT PUNTO 2 (OUTPUT TESTO) -----------------

fprintf('\n============================================================\n');
fprintf('                 PUNTO 2 - FUNZIONE G(s)                    \n');
fprintf('============================================================\n\n');

fprintf('Punto 2: funzione di trasferimento G(s) = ΔT_{out}(s) / ΔPE(s)\n\n');

fprintf('Guadagno statico:\n');
fprintf('  G(0) = %.6f  [°C/W]\n\n', G0);

fprintf('Poli di G(s):\n');
for k = 1:length(p)
    fprintf('  p_%d = %+10.4f  [1/s]\n', k, p(k));
end

if isempty(z)
    fprintf('\nZeri di G(s): nessuno (nessuno zero finito).\n\n');
else
    fprintf('\nZeri di G(s):\n');
    for k = 1:length(z)
        fprintf('  z_%d = %+10.4f  [1/s]\n', k, z(k));
    end
    fprintf('\n');
end

fprintf('Si osserva che i poli di G(s) coincidono con gli autovalori della matrice A.\n');
fprintf('Tutti i poli presentano parte reale negativa; il sistema è quindi BIBO stabile.\n\n');


fprintf('Specifiche dinamiche e parametri derivati:\n');
fprintf('  Sovraelongazione massima S_max       = %.1f  [%%]\n', Sovr_max*100);
fprintf('  Smorzamento equivalente  ξ*          = %.3f  [-]\n', csi_star);
fprintf('  Margine di fase minimo   M_f,des     = %.1f  [deg]\n', Mf_des);
fprintf('  Stima pulsazione di taglio ω_c,min   = %.3f  [rad/s]\n', omega_c_min);
fprintf('  Banda disturbo: [0, %.2f]           [rad/s]\n', omega_dist_max);
fprintf('  Banda rumore : [%.2e, %.2e]         [rad/s]\n', omega_n_min, omega_n_max);
fprintf('------------------------------------------------------------\n\n');


%% ----------------- BODE DI G(s) CON BANDE E VINCOLI -----------------

fig = figure('Name','Punto 2 - Bode di G(s)','NumberTitle','off');

tlo = tiledlayout(2,1);
tlo.TileSpacing = 'compact';
tlo.Padding     = 'compact';

% ----------------- MODULO -----------------
ax1 = nexttile(tlo,1);
hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'XScale','log');
xlim(ax1,[omega_plot_min omega_plot_max]);
ylim(ax1,[mag_min mag_max]);

% Banda disturbo [0, omega_dist_max] (grigio)
hDist_mod = patch(ax1, [omega_plot_min, omega_dist_max, omega_dist_max, omega_plot_min], [mag_max,        mag_max,       mag_min,        mag_min],  [0.8 0.8 0.8], 'FaceAlpha',0.3,'EdgeColor','none');

% Banda rumore [omega_n_min, omega_n_max] (giallo)
hNoise_mod = patch(ax1,  [omega_n_min, omega_n_max, omega_n_max, omega_n_min], [mag_max,     mag_max,     mag_min,     mag_min], [1 1 0], 'FaceAlpha',0.3,'EdgeColor','none');

% Modulo di G(jω)
hG_mod = plot(ax1, wG, MagG_dB, 'LineWidth',1.5, 'Color',[0 0.3 0.8]);

% Linea verticale in ω_c,min
hWc_mod = xline(ax1, omega_c_min, 'r--', '\omega_{c,min}', 'LineWidth',1.2,'LabelOrientation','horizontal');

ylabel(ax1, '|G(j\omega)| [dB]');
title(ax1, 'Punto 2 - Modulo di G(s)');

legend(ax1, [hG_mod, hDist_mod, hNoise_mod, hWc_mod], {'|G(j\omega)|', 'Banda disturbo', 'Banda rumore', '\omega_{c,min}'}, 'Location','southwest', 'Interpreter','tex', 'FontSize',9);

% ----------------- FASE -----------------
ax2 = nexttile(tlo,2);
hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'XScale','log');
xlim(ax2,[omega_plot_min omega_plot_max]);
ylim(ax2,[phi_min phi_max]);

% Limiti locali per patch
yl = ylim(ax2);
phi_low  = yl(1);
phi_high = yl(2);

% Banda disturbo (grigio)
hDist_phase = patch(ax2, [omega_plot_min, omega_dist_max, omega_dist_max, omega_plot_min], [phi_high,       phi_high,       phi_low,       phi_low], [0.8 0.8 0.8], 'FaceAlpha',0.2,'EdgeColor','none');

% Banda rumore (giallo)
hNoise_phase = patch(ax2, [omega_n_min, omega_n_max, omega_n_max, omega_n_min], [phi_high,    phi_high,    phi_low,    phi_low], [1 1 0], 'FaceAlpha',0.2,'EdgeColor','none');

% Zona "proibita" di fase per ω ≥ ω_c,min (margine di fase insufficiente)
% da fase = Mf_des - 180° verso il basso
hZona_phase = patch(ax2, [omega_c_min, omega_n_min, omega_n_min, omega_c_min],[Mf_des-180,  Mf_des-180,  phi_low,     phi_low], [1 0.7 0.7], 'FaceAlpha',0.5,'EdgeColor','none');

% Fase di G(jω) 
hG_phase = plot(ax2, wG, PhaseG_deg, 'LineWidth',1.5, 'Color',[0 0.3 0.8]);

% Linea verticale in ω_c,min
hWc_phase = xline(ax2, omega_c_min, 'r--', '\omega_{c,min}', 'LineWidth',1.2,'LabelOrientation','horizontal');

xlabel(ax2, '\omega [rad/s]');
ylabel(ax2, 'Fase [deg]');
title(ax2, 'Punto 2 - Fase di G(s)');

legend(ax2,  [hG_phase, hDist_phase, hNoise_phase, hZona_phase, hWc_phase], {'\angle G(j\omega)', 'Banda disturbo',  'Banda rumore', 'Zona proibita M_f',  '\omega_{c,min}'}, 'Location','southwest', 'Interpreter','tex','FontSize',9);


%% ======================= PUNTO 3 ==========================
%   PROGETTO DEL REGOLATORE R(s) E VERIFICA SPECIFICHE
% ===========================================================

s = tf('s');    

%% ========================================================================
%  PARTE 1 - CALCOLI (SPECIFICHE, REGOLATORE, L,S,T, RISPOSTE)
% ========================================================================

%% 3.1 - SPECIFICHE 
e_star     = 0.002;    % [°C]   errore a regime massimo
W_max      = 4;        % ampiezza gradino riferimento
D_max      = 3.5;      % ampiezza massima disturbo
Mf_min     = 50;       % [deg]  margine di fase minimo

T_star     = 0.01;     % [s]    tempo di assestamento max (5%)
Sovr_max   = 0.11;     % [-]    overshoot max (11%)

A_d        = 50;       % [dB]   attenuazione disturbo (bassa freq)
omega_dMAX = 0.4;      % [rad/s] limite sup. banda disturbo

A_n        = 60;       % [dB]   attenuazione rumore (alta freq)
omega_n_min = 8e4;     % [rad/s]
omega_n_max = 9e6;     % [rad/s]

%% 3.2 - SPECIFICHE DERIVATE

% Smorzamento da S%
csi_star = abs(log(Sovr_max)) / sqrt(pi^2 + (log(Sovr_max))^2);

% Margine di fase richiesto
Mf = max(csi_star*100, Mf_min);

% Pulsazione di taglio minima (euristica)
omega_c_min = 300 / (Mf * T_star);

%% 3.3 - REGOLATORE STATICO R_s 

% Vincoli:
%   |L(0)| >= (W_max + D_max)/e_star
%   |L(jω_dMAX)| >= 10^(A_d/20)

mu_s_error = (W_max + D_max) / e_star;
mu_s_dist  = 10^(A_d/20);

G_0       = abs( evalfr(G, 0) );
G_om_dMAX = abs( evalfr(G, 1j*omega_dMAX) );

R_s_err = mu_s_error / G_0;
R_s_dis = mu_s_dist  / G_om_dMAX;
R_s     = max(R_s_err, R_s_dis);

% Impianto con solo guadagno statico
G_e = series(R_s, G);

%% 3.4 - REGOLATORE DINAMICO: polo anti-rumore + rete anticipatrice

% 3.4.1 Polo per attenuazione rumore di misura
k_p           = 0.3;     % fattore per posizionare il polo anti-rumore sotto la banda del rumore
omega_p_noise = k_p * omega_n_min;
T1_f          = 1 / omega_p_noise;

R_d_polo = 1 / (1 + T1_f * s);

% 3.4.2 Rete anticipatrice (Scenario B: progetto della rete lead)
omega_c_star = 1.1 * omega_c_min;

% Open-loop senza rete anticipatrice dinamica: L0(s) ≈ G_e(s) = R_s G(s)
G_e_jw = evalfr(G_e, 1j * omega_c_star);
mag_Ge = abs(G_e_jw);
arg_Ge = rad2deg(angle(G_e_jw));   % fase in [−180,180] 

% Margine di fase attuale approssimato:
PM0 = 180 + arg_Ge;                % [deg]

% Normalizza in intervallo [-180,180)
if PM0 > 180
    PM0 = PM0 - 360;
elseif PM0 <= -180
    PM0 = PM0 + 360;
end

margine_sicurezza = 8;             % [deg] margine addizionale per robustezza
PM_des = Mf + margine_sicurezza;   % margine di fase desiderato

% Fase che la rete deve AGGIUNGERE
phi_star = PM_des - PM0;           % [deg]
phi_rad  = phi_star * pi / 180;

% Modulo richiesto della rete a ω_c*
M_star = 1 / mag_Ge;

% Formule Scenario B
tau       = (M_star - cos(phi_rad)) / (omega_c_star * sin(phi_rad));
alpha_tau = (cos(phi_rad) - 1/M_star) / (omega_c_star * sin(phi_rad));

R_d_antic = (1 + tau * s) / (1 + alpha_tau * s);

% Regolatore dinamico completo
R_d = series(R_d_antic, R_d_polo);

% Regolatore completo
R = series(R_s, R_d);

%% 3.5 - FUNZIONI L(s), S(s), T(s)

L = series(R, G);        % anello aperto
S = 1 / (1 + L);         % sensitività
T = L / (1 + L);         % sensitività complementare
F = T;                   % da riferimento w all'uscita

%% 3.6 - CHECK IN FREQUENZA (S, T) E MARGINI

omega_plot_min = 1e-3;
omega_plot_max = 1e8;
w = logspace(-3,8,2000);

[MagL,PhaseL,~] = bode(L,w);
[MagS,PhaseS,~] = bode(S,w);
[MagT,PhaseT,~] = bode(T,w);

MagL_dB = 20*log10(squeeze(MagL));
MagS_dB = 20*log10(squeeze(MagS));
MagT_dB = 20*log10(squeeze(MagT));

PhaseL_deg = squeeze(PhaseL);
PhaseS_deg = squeeze(PhaseS);
PhaseT_deg = squeeze(PhaseT);

idx_d = w <= omega_dMAX;
idx_n = (w >= omega_n_min) & (w <= omega_n_max);

maxS_d = max(MagS_dB(idx_d));   % target: <= -A_d
maxT_n = max(MagT_dB(idx_n));   % target: <= -A_n

[GM, PM, Wgm, Wpm] = margin(L); % margini di stabilità

%% 3.7 - CHECK NEL TEMPO (GRADINO, DISTURBO, RUMORE)

t_sim = 0:1e-4:0.1;
T_sim = t_sim(end);

% Gradino di riferimento
[y_w, t_w] = step(W_max * F, t_sim);
info = stepinfo(W_max * F, 'SettlingTimeThreshold',0.05);

% Disturbo a bassa frequenza (entro [0,0.4])
d_t = 0.2*sin(0.1*1*t_sim) + 0.2*sin(0.1*2*t_sim) + 0.2*sin(0.1*3*t_sim) + 0.2*sin(0.1*4*t_sim);
y_d = lsim(S, d_t, t_sim);

% Rumore di misura in banda alta [8e4, 9e6]
n_t = 0.2*sin(1e5*t_sim) + 0.2*sin(2e5*t_sim) +0.2*sin(3e5*t_sim) + 0.2*sin(4e5*t_sim);
y_n = lsim(-F, n_t, t_sim);   % rumore sul sensore

% Risposta totale
y  = y_w + y_d + y_n;
LV = evalfr(W_max * F, 0);    % valore di regime su riferimento

% Bode "riassuntivo" per G_e e L
w_plot = logspace(-3, 5, 3000);
[MagGe, PhaseGe] = bode(G_e, w_plot);
[MagL2, PhaseL2] = bode(L,   w_plot);

MagGe_dB      = 20*log10(squeeze(MagGe));
MagL2_dB      = 20*log10(squeeze(MagL2));
PhaseGe_deg2  = squeeze(PhaseGe);
PhaseL_deg2   = squeeze(PhaseL2);

idx_wc = find(MagL2_dB <= 0, 1);
wc = w_plot(idx_wc);          % frequenza di crossover

phase_at_Wpm = interp1(w_plot, PhaseL_deg2, Wpm);


%% ========================================================================
%  PARTE 2 - OUTPUT TESTUALE E GRAFICI
% ========================================================================

%% -------------------- OUTPUT TESTUALE PUNTO 3 --------------------------

fprintf('\n============================================================\n');
fprintf('                 PUNTO 3 - PROGETTO DEL REGOLATORE          \n');
fprintf('============================================================\n\n');

fprintf('[3.2] Specifiche derivate:\n');
fprintf('  csi*                = %.4f\n', csi_star);
fprintf('  Mf_des              = %.2f [deg]\n', Mf);
fprintf('  omega_c,min         = %.3f [rad/s]\n\n', omega_c_min);

fprintf('[3.3] Regolatore statico R_s:\n');
fprintf('  R_s (errore a regime) = %.4e\n', R_s_err);
fprintf('  R_s (disturbo)        = %.4e\n', R_s_dis);
fprintf('  R_s scelto            = %.4e\n\n', R_s);

fprintf('[3.4] Rete anticipatrice (Scenario B):\n');
fprintf('  Fase G_e(jω_c*)      = %.2f [deg]\n', arg_Ge);
fprintf('  PM_0 (stima iniziale)  = %.2f [deg]\n', PM0);
fprintf('  PM_des               = %.2f [deg]\n', PM_des);
fprintf('  phi_star (rete)      = %.2f [deg]\n', phi_star);
fprintf('  M_star (|R_d(jω_c*)|)= %.4f\n', M_star);
fprintf('  tau                  = %.4e [s]\n', tau);
fprintf('  alpha_tau            = %.4e [s]\n\n', alpha_tau);


fprintf('[3.5–3.6] Margini di L(s) e attenuazione S, T:\n');
fprintf('  GM (guadagno)        = %.2f [dB]\n', 20*log10(GM));
fprintf('  PM (margine di fase) = %.2f [deg] a ω = %.2f [rad/s]\n', PM, Wpm);
fprintf('  max |S(jω)| in [0, %.1f]      = %.2f [dB] (target <= -%.0f dB)\n', omega_dMAX, maxS_d, A_d);
fprintf('  max |T(jω)| in [%.1e, %.1e] = %.2f [dB] (target <= -%.0f dB)\n\n', omega_n_min, omega_n_max, maxT_n, A_n);

fprintf('[3.7] Prestazioni nel tempo (gradino):\n');
fprintf('  Ts (5%%)             = %.4e [s] (specifica: Ts <= %.4e [s])\n', info.SettlingTime, T_star);
fprintf('  Overshoot            = %.2f [%%] (specifica: S <= %.2f [%%])\n', info.Overshoot, Sovr_max*100);
fprintf('============================================================\n\n');


%% -------------------- GRAFICI PUNTO 3 -------------------------- 

%% 3.6 - Bode di L(s) con margini (margin)
figure('Name','Punto 3 - Bode di L(s)');
margin(L); grid on;
title('Punto 3 - Bode di L(s) = R(s)G(s) e margini di stabilità');


%% 3.7 - Bode L(s), S(s), T(s) (zone + soglie + linea/punto viola)
figure('Name','Punto 3 - Bode L, S, T');
tlo = tiledlayout(2,1);
tlo.TileSpacing = 'compact';
tlo.Padding     = 'compact';

%% --- PALETTE ---
colL = [0   0.25 1];     % blu scuro
colS = [0   0.6  0.2];   % verde elegante
colT = [0.85 0.2  0.2];  % rosso tenue

zoneS = [0.8 1 1];       % disturbo (azzurro )
zoneT = [1 0.9 0.6];     % rumore (arancio )

%% ======================== MODULO ===========================
ax1 = nexttile(tlo,1);
hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'XScale','log');

ylabel(ax1,'Magnitude [dB]');
subtitle_str = sprintf('Gm = %.1f dB (at %.2e rad/s),  Pm = %.1f deg (at %.2f rad/s)', 20*log10(GM), Wgm, PM, Wpm);

title(ax1, {'Bode L(s), S(s), T(s)'; subtitle_str});


ylim(ax1,[-150 100]);
yl = ylim(ax1);

%% --- ZONE COLORATE  ---
hZoneS = patch(ax1, [w(1) omega_dMAX omega_dMAX w(1)],  [yl(1) yl(1) -A_d -A_d],  zoneS, 'FaceAlpha',0.20, 'EdgeColor','none', 'DisplayName','Zona disturbo');
hZoneT = patch(ax1, [omega_n_min omega_n_max omega_n_max omega_n_min],  [yl(1) yl(1) -A_n -A_n],  zoneT, 'FaceAlpha',0.20, 'EdgeColor','none',  'DisplayName','Zona rumore');

%% --- Curve L, S, T ---
hL = plot(ax1, w, MagL_dB, 'Color',colL, 'LineWidth',1.6, 'DisplayName','L(s)');
hS = plot(ax1, w, MagS_dB, 'Color',colS, 'LineWidth',1.6, 'DisplayName','S(s)');
hT = plot(ax1, w, MagT_dB, 'Color',colT, 'LineWidth',1.6, 'DisplayName','T(s)');

%% --- Linee soglia ---
hSdB = yline(ax1, -A_d, 'Color',[0 0.4 1], 'LineStyle','--', 'LineWidth',1.3, 'DisplayName', sprintf('Soglia S = -%d dB',A_d));
hTdB = yline(ax1, -A_n, 'Color',[1 0.6 0], 'LineStyle','--', 'LineWidth',1.3, 'DisplayName', sprintf('Soglia T = -%d dB',A_n));

%% --- Legenda  ---
legend(ax1, [hL hS hT hZoneS hZoneT hSdB hTdB], {'L(s)','S(s)','T(s)', 'Zona disturbo','Zona rumore', sprintf('Soglia S = -%d dB',A_d), sprintf('Soglia T = -%d dB',A_n)}, 'Location','best', 'FontSize',9);

%% ========================= FASE ============================
ax2 = nexttile(tlo,2);
hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'XScale','log');

xlabel(ax2,'Frequency [rad/s]');
ylabel(ax2,'Phase [deg]');
title(ax2,'Fase di L(s), S(s), T(s)');

ylim(ax2,[-270 120]);
yticks(ax2, -270:45:120);

%% Curve fase
hL_p = plot(ax2, w, PhaseL_deg, 'Color',colL, 'LineWidth',1.6, 'DisplayName','L(s)');
hS_p = plot(ax2, w, PhaseS_deg, 'Color',colS, 'LineWidth',1.6, 'DisplayName','S(s)');
hT_p = plot(ax2, w, PhaseT_deg, 'Color',colT, 'LineWidth',1.6, 'DisplayName','T(s)');

%% Linea -180°
yline(ax2, -180, 'Color',[0.4 0.4 0.4], 'LineStyle','--', 'LineWidth',1.2);

%% Linea viola + punto viola (margine di fase Wpm)
if exist('Wpm','var') && ~isempty(Wpm) && ~isnan(Wpm)
    hWpm = xline(ax2, Wpm, 'Color',[0.5 0 0.7], 'LineStyle','--', 'LineWidth',1.6,'DisplayName','\omega_{pm}');
    phase_L_Wpm = interp1(w, PhaseL_deg, Wpm, 'linear','extrap');
    hMarker = plot(ax2, Wpm, phase_L_Wpm, 'o', 'Color',[0.5 0 0.7], 'MarkerSize',7, 'LineWidth',1.5, 'DisplayName','Punto M_f');
end

%% Legenda fase 
legend(ax2, [hL_p hS_p hT_p hWpm hMarker], {'L(s)','S(s)','T(s)', '\omega_{pm}', 'Punto M_f'}, 'Location','best', 'FontSize',9);

%% 3.7 bis - Bode riassuntivo L(s) e G_e(s)

figure('Name','Punto 3 - Bode L(s) Ge(s)');
tlo = tiledlayout(2,1);
tlo.TileSpacing = 'compact';
tlo.Padding     = 'compact';



%% ===================== MODULO =====================
ax1 = nexttile(tlo,1);
hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'XScale','log');

ylabel(ax1,'Ampiezza [dB]');
main_title = 'Bode L(s) e G_e(s)';
subtitle_str = sprintf('Gm = %.1f dB (at %.2e rad/s),  Pm = %.1f deg (at %.2f rad/s)', 20*log10(GM), Wgm, PM, Wpm);

title(ax1, {main_title; subtitle_str});

% colori
colGe = [0   0.25 1];    % blu
colL  = [0.85 0.2  0.2]; % rosso

% --- CURVE ---
hGe = plot(ax1, w_plot, MagGe_dB, 'Color',colGe,'LineWidth',1.8,'DisplayName','G_e(s)');
hL  = plot(ax1, w_plot, MagL2_dB, 'Color',colL, 'LineWidth',1.8,'DisplayName','L(s)');

% linea 0 dB
yline(ax1, 0, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

% scala verticale
ylim(ax1,[-150 100]);
yl = ylim(ax1);

% --- ZONA ω_c ragionevole: da ω_c,min a circa 10·ω_c,min ---
y_top    = min( 20,  yl(2) );
y_bottom = max(-150, yl(1) );

x_left  = omega_c_min;
x_right = min(10*omega_c_min, w_plot(end));   % non più fino a w_plot(end)

hWc_box = patch(ax1, [x_left  x_right  x_right  x_left],[y_top   y_top    y_bottom y_bottom],[1 0.6 0], 'FaceAlpha',0.15,'EdgeColor','none');

% --- Punto di crossover ω_c (|L| = 1 -> 0 dB) ---
L_wc = MagL2_dB(idx_wc);          % MagL2_dB calcolato su w_plot
h_wc_point = plot(ax1, wc, L_wc, 'o', 'MarkerSize',7, 'LineWidth',1.6, 'Color',[1 0.4 0]);

% legenda modulo
legend(ax1, [hGe, hL, hWc_box, h_wc_point], { 'G_e(s)', 'L(s)','\omega \le \omega_{c,min}', '\omega_c: |L| = 1' }, 'Location','eastoutside', 'Interpreter','tex', 'FontSize',9);

%% ===================== FASE =====================
ax2 = nexttile(tlo,2);
hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'XScale','log');
xlabel(ax2,'Frequency [rad/s]');
ylabel(ax2,'Fase [deg]');
title(ax2,'Fase di G_e(s) e L(s)');

% curve di fase
hGe_p = plot(ax2, w_plot, PhaseGe_deg2,'Color',colGe,'LineWidth',1.8);
hL_p  = plot(ax2, w_plot, PhaseL_deg2, 'Color',colL, 'LineWidth',1.8);

% linea -180°
yline(ax2,-180,'--','Color',[.4 .4 .4],'LineWidth',1.2);

% ====== RETTANGOLO: margine di fase EFFETTIVO ======
phi_low  = -180;          % base
phi_high = phase_at_Wpm;  % = -180 + PM

x_left  = omega_c_min;                 % da w_c,min
x_right = min(10*omega_c_min, w_plot(end));

hMf_rect = patch(ax2, [x_left   x_right  x_right  x_left], [phi_low  phi_low  phi_high phi_high], [1 0 0], 'FaceAlpha',0.15,'EdgeColor','none');

% linea soglia per M_f richiesto
phi_req = -180 + Mf;
hMf_line = yline(ax2, phi_req, 'Color',[1 0.6 0],'LineStyle','--', 'LineWidth',1.4);

% ====== PALLINA VIOLA: margine di fase a wpm ======
hMf_point = plot(ax2, Wpm, phase_at_Wpm, 'o', 'Color',[0.7 0 0.9], 'MarkerSize',9, 'LineWidth',1.8);

% ====== PUNTO BLU: margine di guadagno (ωgm) ======
[~, idx_180] = min(abs(PhaseL_deg2 + 180));
w_180   = w_plot(idx_180);
phi_180 = PhaseL_deg2(idx_180);
h_phase180 = plot(ax2, w_180, phi_180, 'o', 'MarkerSize',7, 'LineWidth',1.6, 'Color',[0 0.5 1]);

% legenda fase
legend(ax2, [hMf_rect, hMf_line, hMf_point, h_phase180],{ sprintf('PM = %.1f°', PM), sprintf('Soglia M_f = %.1f°', Mf), '\omega_{pm}', '\omega_{gm} (fase \approx -180°)' }, 'Location','eastoutside', 'Interpreter','tex', 'FontSize',9);

