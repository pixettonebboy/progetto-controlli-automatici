clear ; clc ; close all;

hR     = 50;       
AR     = 0.07;     
cR     = 840.8;    
cA     = 1010;     
mR     = 2.542;   
mA     = 0.1041;   
m_dotA = 0.2;      
Tin    = 28;       
kappa  = 3e-3;    
aR   = hR * AR;        
den1 = mR * cR;       
den2 = mA * cA;       

%% PUNTO 1: MOD>ELLO DI STATO E LINEARIZZAZIONE
% Stato:    x = [TR; Tout]
% Ingresso: u = PE
% Uscita:   y = Tout

% PUNTO DI EQUILIBRIO 
TR_e   = 175;          % °C
Tout_e = 80;           % °C
xe     = [TR_e; Tout_e];

% Fattore KTR = 1 + kappa*TR 
KTR = 1 + kappa * TR_e;

% Calcolo dell'ingresso di equilibrio u_e
% In equilibrio: dTR/dt = 0, dTout/dt = 0.
%
% Dalla prima equazione:
%   0 = (aR/den1)*(Tout_e - TR_e) + u_e / (KTR * den1)
%
% => u_e = - aR * (Tout_e - TR_e) * KTR

u_e = - aR * (Tout_e - TR_e) * KTR;   
ue  = u_e;                            

fprintf('=== PUNTO 1: PUNTO DI EQUILIBRIO ===\n');
fprintf('  TR_e   = %.4f °C\n', TR_e);
fprintf('  Tout_e = %.4f °C\n', Tout_e);
fprintf('  u_e    = %.4f W\n\n', u_e);

% equazioni dinamiche
% f1(TR, Tout, u) = dTR/dt   = (aR/den1)*(Tout - TR) + u / ((1 + kappa*TR)*den1)
% f2(TR, Tout)    = dTout/dt = (m_dotA/mA)*(Tin - Tout) + (aR/den2)*(TR - Tout)
% Linearizzazione: A = df/dx |_(TR_e, Tout_e, u_e), B = df/du |_(TR_e,...)
%   A11 = ∂f1/∂TR   A12 = ∂f1/∂Tout
%   A21 = ∂f2/∂TR   A22 = ∂f2/∂Tout

% Derivate parziali di f1 (dTR/dt)
% f1 = (aR/den1)*(Tout - TR) + u / ( (1 + kappa*TR) * den1 )
%
% ∂f1/∂TR   = -aR/den1  + ∂/∂TR [ u / ( (1 + kappa*TR) * den1 ) ]
%           = -aR/den1  - (u * kappa) / ( den1 * (1 + kappa*TR)^2 )
% Valutata in (TR_e, u_e):
A11 = - aR/den1 - (u_e * kappa) / (den1 * KTR^2);

% ∂f1/∂Tout = (aR/den1)
A12 =  aR/den1;

% Derivate parziali di f2 (dTout/dt)
% f2 = (m_dotA/mA)*(Tin - Tout) + (aR/den2)*(TR - Tout)
%
% ∂f2/∂TR   = (aR/den2)
A21 = aR/den2;

% ∂f2/∂Tout = -(m_dotA/mA) - (aR/den2)
A22 = - m_dotA/mA - aR/den2;

% Matrice A del sistema linearizzato
A = [A11  A12;
     A21  A22];

%  Matrice B (derivate rispetto all'ingresso u = PE)
% f1 dipende da u, f2 NON dipende da u
%
% ∂f1/∂u = 1 / ( (1 + kappa*TR) * den1 ) valutata in TR_e:
B1 = 1 / (den1 * KTR);

% ∂f2/∂u = 0
B2 = 0;

B = [B1;
     B2];

% Matrici C e D (uscita y = Tout)
% y = [0 1] * [TR; Tout]  =>  C = [0 1], D = 0
C = [0 1];
D = 0;
lambda = eig(A);
fprintf('\n=== Verifica stabilità (Punto 1) ===\n');
disp('Autovalori della matrice A:');
disp(lambda);

if all(real(lambda) < 0)
    disp('Il sistema linearizzato è asintoticamente stabile (Re(λ) < 0).');
else
    disp('Il sistema NON è stabile (almeno un autovalore ha Re(λ) ≥ 0)');
end


% Stampa matrici di stato
disp('PUNTO 1: MATRICI DEL SISTEMA LINEARIZZATO');
disp('Matrice A:');
disp(A);
disp('Matrice B:');
disp(B);
disp('Matrice C:');
disp(C);
disp('Matrice D:');
disp(D);

%  autovalori di A
lambda = eig(A);
disp('Autovalori della matrice A (stabilita'' del modello lineare):');
disp(lambda);
fprintf('\n');
