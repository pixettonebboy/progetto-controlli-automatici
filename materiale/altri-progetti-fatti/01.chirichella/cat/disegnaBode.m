
function disegnaBode(G, titolo, ZP_d, ZP_n, ZP_OmegaC, ZP_Mf)
    
    % Funzione per creare un diagramma di Bode con specifiche e annotazioni
    omega_plot_min = 1e-3;
    omega_plot_max = 1e8;

    % Creazione figura Bode
    figure('Name', titolo, 'NumberTitle', 'off');
    hold on; grid on; zoom on;
    
    % Specifiche su d (attenuazione del disturbo sull'uscita)
    patch(ZP_d{1}, ZP_d{2}, 'r', 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'DisplayName', 'Attenuazione Disturbo d');
    
    % Specifiche su n (attenuazione del disturbo di misura)
    patch(ZP_n{1}, ZP_n{2}, 'g', 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'DisplayName', 'Attenuazione Disturbo n');
    
    % Specifiche su omega_c,min (frequenza critica minima)
    patch(ZP_OmegaC{1}, ZP_OmegaC{2}, 'b', 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'DisplayName', '\omega_{c,min}');
    
    % Diagramma di Bode per L con margini di stabilità
    margin(G, {omega_plot_min, omega_plot_max});
    
    % Specifiche sul margine di fase
    patch(ZP_Mf{1}, ZP_Mf{2}, 'm', 'FaceAlpha', 0.2, 'EdgeAlpha', 0, 'DisplayName', 'Margine di fase proibito');
    
    % Legenda colori
    legend('show', 'Location', 'southwest');
    
    % Titoli e annotazioni
    title(titolo);
    xlabel('Frequenza [rad/s]');
    ylabel('Ampiezza [dB] e Fase [°]');
    
    hold off;
end
