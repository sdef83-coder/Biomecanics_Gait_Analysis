function sparc = compute_SPARC(velocity, fs, varargin)
% compute_SPARC - Calcule le SPARC avec option de visualisation du spectre
%
% Syntaxe:
%   sparc = compute_SPARC(velocity, fs)
%   sparc = compute_SPARC(velocity, fs, 'plot', true)
%   sparc = compute_SPARC(velocity, fs, 'plot', true, 'title', 'Mon titre')
%
% Paramètres:
%   velocity : vecteur de vitesse
%   fs       : fréquence d'échantillonnage (Hz)
%   'plot'   : true/false pour afficher le spectre (défaut: false)
%   'title'  : titre personnalisé pour la figure

    % Parsing des arguments optionnels
    p = inputParser;
    addParameter(p, 'plot', false, @islogical);
    addParameter(p, 'title', '', @ischar);
    parse(p, varargin{:});

    doPlot = p.Results.plot;
    plotTitle = p.Results.title;

    % Vérifications initiales
    if length(velocity) < 20 || all(isnan(velocity))
        sparc = NaN;
        return;
    end

    velocity = velocity(:);

    % Etapes basées sur Balasubramanian et al. 2012

    % ÉTAPE 1: FFT avec zero-padding (padlevel = 4)
    N = length(velocity);
    K = 2^(nextpow2(N) + 4);  
    V = abs(fft(velocity, K));
    V = V(1:K/2+1);

    % ÉTAPE 2: Normalisation par DC
    V0 = max(V(1), 1e-10);
    Vn = V / V0; 
    %Vn = V / max(V); % test normalisation par max du spectre (pas recommandé)

    % ÉTAPE 3: Arc length sur 0-6 Hz
    freqs = (0:K/2)' * (fs / K);
    f_c = 6;
    omega = 2 * pi * freqs;
    omega_c = 2 * pi * f_c;  % 6 Hz
    idx = find(omega <= omega_c);

    if length(idx) < 2
        sparc = NaN;
        return;
    end

    omega_seg = omega(idx);
    V_seg = Vn(idx);
    d_omega = diff(omega_seg);
    d_V = diff(V_seg);
    arc = sum(sqrt((d_omega / omega_c).^2 + d_V.^2));
    sparc = -arc;

    % VISUALISATION
    if doPlot
        figure('Name', 'Analyse spectrale SPARC', 'Color', 'w');

        % Subplot 1: Signal temporel
        subplot(3,1,1);
        t = (0:length(velocity)-1) / fs;
        plot(t, velocity, 'b-', 'LineWidth', 1.5);
        grid on; box on;
        xlabel('Temps (s)');
        ylabel('Vitesse (m/s)');
        title('Signal temporel');

        % Subplot 2: Spectre complet (0 à Nyquist)
        subplot(3,1,2);
        plot(freqs, Vn, 'k-', 'LineWidth', 1);
        hold on;
        xline(6, 'r--', 'LineWidth', 2, 'Label', '6 Hz (cutoff)');
        hold off;
        grid on; box on;
        xlabel('Fréquence (Hz)');
        ylabel('Magnitude normalisée');
        title('Spectre fréquentiel complet');
        xlim([0 min(50, fs/2)]);

        % Subplot 3: Spectre dans la bande SPARC (0-6 Hz) avec arc-length
        subplot(3,1,3);
        freq_seg = omega_seg / (2*pi);
        plot(freq_seg, V_seg, 'b-', 'LineWidth', 2);
        hold on;

        % Visualisation de l'arc-length
        plot(freq_seg, V_seg, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');

        hold off;
        grid on; box on;
        xlabel('Fréquence (Hz)');
        ylabel('Magnitude normalisée');
        title(sprintf('Bande SPARC (0-6 Hz) | SPARC = %.3f | Arc-length = %.3f', sparc, arc));
        xlim([0 6]);

        % Titre global
        if ~isempty(plotTitle)
            sgtitle(plotTitle, 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
end