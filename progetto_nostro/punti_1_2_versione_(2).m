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
%       LINEARIZZAZIONE: CALCOLO NUMERICO DI A_e E B_e (stile riferimento)
% -------------------------------------------------------------------------

delta = 1e-6;   % incremento per derivate finte

% Matrice A_e: derivate rispetto a [TR, Tout]
A_e = zeros(2,2);

% df1/dTR
A_e(1,1) = ( f1(x1_e + delta, x2_e, u_e) - f1(x1_e, x2_e, u_e) ) / delta;
% df1/dTout
A_e(1,2) = ( f1(x1_e, x2_e + delta, u_e) - f1(x1_e, x2_e, u_e) ) / delta;

% df2/dTR
A_e(2,1) = ( f2(x1_e + delta, x2_e, u_e) - f2(x1_e, x2_e, u_e) ) / delta;
% df2/dTout
A_e(2,2) = ( f2(x1_e, x2_e + delta, u_e) - f2(x1_e, x2_e, u_e) ) / delta;

% Matrice B_e: derivate rispetto a u
B_e = zeros(2,1);

% df1/du
B_e(1,1) = ( f1(x1_e, x2_e, u_e + delta) - f1(x1_e, x2_e, u_e) ) / delta;
% df2/du
B_e(2,1) = ( f2(x1_e, x2_e, u_e + delta) - f2(x1_e, x2_e, u_e) ) / delta;

% Uscita: y = Tout
C = [0 1];
D = 0;

% Verifica stabilità
lambda = eig(A_e);
fprintf('Autovalori della matrice A_e:\n');
disp(lambda);
if all(real(lambda) < 0)
    disp('Il sistema linearizzato è asintoticamente stabile (Re(λ) < 0).');
else
    disp('Il sistema linearizzato NON è asintoticamente stabile.');
end

% Costruiamo il modello linearizzato
sys_ss = ss(A_e, B_e, C, D);


%% ------------------------------------------------------------------------
%                         PUNTO 2: FUNZIONE G(s)
% -------------------------------------------------------------------------

% Funzione di trasferimento G(s) = ΔTout(s) / ΔPE(s)
G = tf(sys_ss);

fprintf('\n-PUNTO 2: Funzione di trasferimento G(s) -\n');
zpk(G)    % forma zero-poli-guadagno
G0 = dcgain(G);
fprintf('Guadagno statico G(0) = %.6f [°C/W]\n', G0);
