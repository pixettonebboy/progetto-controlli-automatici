%% Funzione di trasferimento del sistema linearizzato
% Costruzione del modello
s  = tf('s');
GG = C*inv(s*eye(2) - A)*B + D;

%% Diagramma di bode della funzione di trasferimento
figure('Name', 'Funzione di trasferimento');
hold on;
bode(GG,{omega_plot_min,omega_plot_max});
title_str = "Funzione di trasferimento" + newline + "Diagaramma di Bode";
title(title_str);
grid on, zoom on;
hold off;
