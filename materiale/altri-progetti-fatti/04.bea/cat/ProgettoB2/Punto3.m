%% Sintesi regolatore statico
% Valore minimo prescritto per L(0)
mu_s_error  = (DD+WW)/e_star;       
mu_s_dist   = 10^(A_d/20);          

% Calcolo valori in 0 e Omega MAX
G_0         = abs(evalfr(GG, 0));
G_om_d_MAX  = abs(evalfr(GG, 1j*omega_d_MAX));

% Definizione regolatore statico
v1 = mu_s_error/G_0;
v2 = mu_s_dist/G_om_d_MAX;
RR_s        = max(mu_s_error/G_0, mu_s_dist/G_om_d_MAX); 

fprintf('v1: %f\n', v1);
fprintf('v2: %f\n', v2);
fprintf('Regolatore Statico: %f\n', RR_s);

% Sistema esteso
GG_e        = RR_s*GG;

%% Diagrammi di Bode del sistema esteso con specifiche
figure('Name', 'Plot sistema esteso');
hold on;

% Calcolo specifiche S% => Margine di fase
xi_star     = abs(log(S_star/100))/sqrt(pi^2 + log(S_star/100)^2);
Mf          = max(xi_star*100, Mf_esp);

% Specifiche su d
Bnd_d_x     = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y     = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x     = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y     = [-A_n; -A_n; 100; 100];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione critica)
omega_Ta_min = omega_plot_min; % lower bound per il plot
omega_Ta_max = 460/(Mf*T_star); % omega_c >= 460/(Mf*T^*)
Bnd_Ta_x     = [omega_Ta_min; omega_Ta_max; omega_Ta_max; omega_Ta_min];
Bnd_Ta_y     = [0; 0; -150; -150];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);

% Legenda colori
%Legend_mag  = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
%legend(Legend_mag);

title_str = "Diagaramma di Bode";


% Plot Bode con margini di stabilità
margin(GG_e,{omega_plot_min,omega_plot_max});
title(title_str);
grid on; zoom on;

% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_max;
omega_c_max = omega_n_min;

phi_up      = Mf - 180;
phi_low     = -270; % lower bound per il plot

Bnd_Mf_x    = [omega_c_min; omega_c_max; omega_c_max; omega_c_min];
Bnd_Mf_y    = [phi_up; phi_up; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

%Legenda
LegendaModulo = ["A_d"; "GG_e"];
%legend(LegendaModulo);

% Legenda colori
Legend_mag  = ["\omega_{c,min}"; "A_n"; "A_d"; "G(j\omega)"];
legend(Legend_mag);

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

hold off;

%% Design del regolatore dinamico
% Loopshaping per definire omega_c_star e Mf_star

Mf_star             = Mf+10;
omega_c_star        = 1900;

mag_omega_c_star_dB = abs(evalfr(GG_e, 1j*omega_c_star));
arg_omega_c_star    = rad2deg(angle(evalfr(GG_e, 1j*omega_c_star)));

M_star      = 1/mag_omega_c_star_dB;        % M_star = 1/GG_e(omega_c_star)
phi_star    = Mf_star - 180 - arg_omega_c_star;

tau         = (M_star - cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha_tau   = (cos(phi_star*pi/180) - 1/M_star)/(omega_c_star*sin(phi_star*pi/180));
alpha       = alpha_tau / tau;

fprintf('tau: %f\n', tau);
fprintf('alpha: %f\n', alpha);

if min(tau,alpha) < 0
    fprintf('Errore: parametri rete anticipatrice negativi\n');
    return;
end


RR_d = (1 + tau*s)/(1 + alpha*tau*s); 


RR = RR_s*RR_d;

LL = RR*GG; % funzione di anello


%% Diagrammi di Bode con regolatore dinamico

figure('Name', 'Plot con rete anticipatrice');
hold on;


% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL,{omega_plot_min,omega_plot_max});
title(title_str);
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_arg);
hold off;
