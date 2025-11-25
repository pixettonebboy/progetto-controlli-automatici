%% Definizione delle variabili
% u     ->  Coppia in ingresso
% x1    ->  theta (posizione del sistema)
% x2    ->  omega (velocità del sistema)
% y     ->  theta (posizione del sistema)
syms u x1 x2 y;

% Definizioni funzioni di stato
f = cos(alpha)/1-(sin(alpha)*cos(x1))^2;
f1 = x2;
f2 = (u*cos(alpha))/(J*(1-(sin(alpha)*cos(x1))^2)) - (beta*x2)/J - (k*x1)/J;

% Definizione funzione di uscita
h = x1;

%% Calcolo equilibrio
x1e = deg2rad(140);     % Da testo
x2e = 0;                % f1(xe, ue) = 0
ye  = deg2rad(140);     
ue  = k*x1e / (cos(alpha)/(1-(sin(alpha)*cos(x1e))^2));      

%% Linearizzazione del sistema
% Calcolo delle Jacobiane
A_raw = jacobian([f1, f2], [x1, x2]);
B_raw = jacobian([f1, f2], u);
C_raw = jacobian(h, [x1, x2]);
D_raw = jacobian(h, u);

% Valori nel punto di equilibrio
A = double(subs(A_raw, [x1, x2, u], [x1e, x2e, ue]));
B = double(subs(B_raw, [x1, x2, u], [x1e, x2e, ue]));
C = double(subs(C_raw, [x1, x2, u], [x1e, x2e, ue]));
D = double(subs(D_raw, [x1, x2, u], [x1e, x2e, ue]));

disp(A);
disp(B);
disp(C);
disp(D);