clear all; close all; clc;

%% ------------------------------------------------------------------------
%                           FASE DI ANALISI
% -------------------------------------------------------------------------

%% PARAMETRI DEL SISTEMA 
hR = 50;        % [W/(m^2 C°)]
AR = 0.07;      % [m^2]
cR = 840.8;     % [J/(kg C°)]
cA = 1010;      % [J/(kg C°)]
mR = 2.542;     % [kg]
mA = 0.1041;    % [kg]
m_dotA = 0.2;       % [kg/s]
Tin = 28;        % [C°]
kappa = 3e-3;      % [1/C°]

aR = hR * AR;            % coeff. convezione complessivo [W/C°]
den1 = mR * cR;            % mR*cR
den2 = mA * cA;            % mA*cA


%% DEFINIZIONE FUNZIONI DI STATO
% x = [TR, Tout],  u = PE (potenza elettrica)
% risolto per x1' = dTR/dt e x2' = dTout/dt

% Equazione (1):
% mR*cR * dTR/dt = hR*AR*(Tout - TR) + PE/(1 + kappa*TR)
f1 = @(TR, Tout, u)( aR*(Tout - TR) + u./(1 + kappa*TR) ) / den1;

% Equazione (2):
% mA*cA * dTout/dt = m_dotA*cA*(Tin - Tout) + hR*AR*(TR - Tout)
f2 = @(TR, Tout, u)( m_dotA*cA*(Tin - Tout) + aR*(TR - Tout) ) / den2;


%% ------------------------------------------------------------------------
%                        EQUILIBRIO DEL SISTEMA (Punto 1)
% -------------------------------------------------------------------------

% Valori di equilibrio (da tabella)
TR_e = 175;          % [C°]
Tout_e = 80;           % [C°]
x1_e = TR_e;
x2_e = Tout_e;
x_e  = [ x1_e ; x2_e ];

% Calcolo ingresso di equilibrio u_e = PE_e
% Dalla wquazione (1) in equilibrio: 0 = aR*(Tout_e - TR_e) + u_e/(1 + kappa*TR_e)
KTR_e = 1 + kappa * TR_e;
u_e  = - aR * (Tout_e - TR_e) * KTR_e;   % [W]

fprintf('--- PUNTO 1: EQUILIBRIO ---\n');
fprintf('TR_e   = %.4f °C\n', TR_e);
fprintf('Tout_e = %.4f °C\n', Tout_e);
fprintf('u_e    = %.4f W\n\n', u_e);


%% ------------------------------------------------------------------------
%       LINEARIZZAZIONE: CALCOLO NUMERICO DI A_e E B_e 
% -------------------------------------------------------------------------

delta = 1e-6;   % incremento per derivate finte

% Matrice A: derivate rispetto a [TR, Tout]
A = zeros(2,2);

% df1/dTR
A(1,1) = ( f1(x1_e + delta, x2_e, u_e) - f1(x1_e, x2_e, u_e) ) / delta;
% df1/dTout
A(1,2) = ( f1(x1_e, x2_e + delta, u_e) - f1(x1_e, x2_e, u_e) ) / delta;

% df2/dTR
A(2,1) = ( f2(x1_e + delta, x2_e, u_e) - f2(x1_e, x2_e, u_e) ) / delta;
% df2/dTout
A(2,2) = ( f2(x1_e, x2_e + delta, u_e) - f2(x1_e, x2_e, u_e) ) / delta;

% Matrice B: derivate rispetto a u
B = zeros(2,1);

% df1/du
B(1,1) = ( f1(x1_e, x2_e, u_e + delta) - f1(x1_e, x2_e, u_e) ) / delta;
% df2/du
B(2,1) = ( f2(x1_e, x2_e, u_e + delta) - f2(x1_e, x2_e, u_e) ) / delta;

% Uscita: y = Tout
C = [0 1];
D = 0;

% Verifica stabilità
lambda = eig(A);
fprintf('Autovalori della matrice A:\n');
disp(lambda);
if all(real(lambda) < 0)
    disp('Il sistema linearizzato è asintoticamente stabile (Re(λ) < 0).');
else
    disp('Il sistema linearizzato NON è asintoticamente stabile.');
end

%  modello linearizzato
sys_ss = ss(A, B, C, D);


%% ------------------------------------------------------------------------
%                         PUNTO 2: FUNZIONE G(s)
% -------------------------------------------------------------------------

%% ======================= PUNTO 2 ==========================
% Funzione di trasferimento G(s) = ΔTout / ΔPE
% Bode di G(s) con:
%  - banda disturbo (grigio)
%  - banda rumore (giallo)
%  - linea in ω_c,min
%  - zona proibita per Mf
% -------------------------------------------------------------------------

sys_ss = ss(A, B, C, D);
G = tf(sys_ss);

fprintf('--- PUNTO 2: G(s) = ΔTout(s) / ΔPE(s) ---\n');
zpk(G)
G0 = dcgain(G);
fprintf('Guadagno statico G(0) = %.6f [°C/W]\n\n', G0);

%  SPECIFICHE PER LE ZONE
Mf_min   = 50;       % [deg]
T_star   = 0.01;     % [s]
Sovr_max = 0.11;     % 11%

A_d      = 50;       % [dB] attenuazione disturbo in [0, 0.4] rad/s
A_n      = 60;       % [dB] attenuazione rumore in [8e4, 9e6] rad/s

omega_dist_max = 0.4;     % banda disturbo
omega_n_min    = 8e4;     % inizio banda rumore
omega_n_max    = 9e6;     % fine banda rumore

% Smorzamento e Mf desiderato
csi_star = abs(log(Sovr_max)) / sqrt(pi^2 + (log(Sovr_max))^2);
Mf       = max(csi_star*100, Mf_min);

% ω_c,min euristica
omega_c_min = 300 / (Mf * T_star);

fprintf('PUNTO 2: Mf_des = %.2f gradi, omega_c,min = %.3f rad/s\n\n', Mf, omega_c_min);

% BODE  PER MODULO E FASE 
omega_plot_min = 1e-3;
omega_plot_max = 1e7;

[MagG, PhaseG, wG] = bode(G, {omega_plot_min, omega_plot_max});
MagG       = squeeze(MagG);
PhaseG_deg = squeeze(PhaseG);          % [deg]
MagG_dB    = 20*log10(MagG);

%% BODE COMPLETO 
fig = figure('Name','Punto 2 - Bode di G(s) con bande e vincoli');
tlo = tiledlayout(2,1);  

% ----------------- MODULO -----------------
ax1 = nexttile(tlo,1);
hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
set(ax1,'XScale','log');

% banda disturbo [0, omega_dist_max]
patch(ax1, [omega_plot_min, omega_dist_max, omega_dist_max, omega_plot_min], [80, 80, -140, -140], [0.8 0.8 0.8], 'FaceAlpha',0.3,'EdgeColor','none');

% banda rumore [omega_n_min, omega_n_max]
patch(ax1, [omega_n_min, omega_n_max, omega_n_max, omega_n_min], [80, 80, -140, -140], [1 1 0], 'FaceAlpha',0.3,'EdgeColor','none');

plot(ax1, wG, MagG_dB, 'LineWidth',1.5);
xline(ax1, omega_c_min, 'r--', '\omega_{c,min}', 'LineWidth',1.2,'LabelOrientation','horizontal');

ylabel(ax1, '|G(j\omega)| [dB]');
title(ax1, 'Punto 2 - Bode di G(s) con bande disturbo/rumore e \omega_{c,min}');

% ----------------- FASE -----------------
ax2 = nexttile(tlo,2);
hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
set(ax2,'XScale','log');

plot(ax2, wG, PhaseG_deg, 'LineWidth',1.5);

yl = ylim(ax2);
phi_low = yl(1);

% banda disturbo
patch(ax2, [omega_plot_min, omega_dist_max, omega_dist_max, omega_plot_min], [yl(2), yl(2), yl(1), yl(1)],[0.8 0.8 0.8],'FaceAlpha',0.2,'EdgeColor','none');

% banda rumore
patch(ax2, [omega_n_min, omega_n_max, omega_n_max, omega_n_min], [yl(2), yl(2), yl(1), yl(1)],[1 1 0],'FaceAlpha',0.2,'EdgeColor','none');

% zona proibita per Mf (fase < Mf-180 a destra di ω_c,min)
patch(ax2, [omega_c_min, omega_n_min, omega_n_min, omega_c_min],[Mf-180, Mf-180, phi_low, phi_low],[1 0.7 0.7], 'FaceAlpha',0.5,'EdgeColor','none');

xline(ax2, omega_c_min, 'r--', '\omega_{c,min}','LineWidth',1.2,'LabelOrientation','horizontal');

xlabel(ax2, '\omega [rad/s]');
ylabel(ax2, 'Fase [deg]');
title(ax2, 'Fase di G(s) con bande e zona proibita Mf');

tlo.TileSpacing = 'compact';
tlo.Padding     = 'compact';


