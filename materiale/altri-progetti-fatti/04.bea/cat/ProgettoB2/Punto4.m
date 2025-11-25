%% Test sul sistema sistema linearizzato

% Funzione di sensitività
SS = 1/(1+LL);

% Funzione di sensitività complementare
FF = LL/(1+LL);

%% Risposta al gradino
figure('Name', 'Simulazione risposta a gradino');

T_simulation = 2*T_star;
%T_simulation = 2*T_star;
[y_step,t_step] = step(WW*FF, T_simulation);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

LV = evalfr(WW*FF,0);

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0], [LV*(1+S_star/100), LV*(1+S_star/100),LV*2,LV*2],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);

% vincolo tempo di assestamento al 5%
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1-0.05),LV*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_star,T_simulation,T_simulation,T_star],[LV*(1+0.05),LV*(1+0.05),LV*2,LV*2],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

ylim([0,LV*2]);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

