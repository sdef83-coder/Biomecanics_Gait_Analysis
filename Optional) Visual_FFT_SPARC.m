%% ================================================================
%  VISUALISATION SPECTRALE SPARC - STYLE ARTICLE BECK ET AL. 2018
%  - Spectres normalisés pour le SPARC (Magnitude COM)
%  - Superposition des 3 conditions (Plat/Medium/High) pour l'essai 1
%  - 2 figures : (1) Spectre non-normalisé, (2) Spectre normalisé
% ================================================================
clear; clc; close all;

%% === CHEMINS ET PARAMÈTRES ===
addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\Scripts\btk'));    
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));

sujet_id  = 'CTL_01';
surfaces  = {'Plat', 'Medium', 'High'};
essai_cible = 1;  % On prend l'essai 1 pour la comparaison
base_dir  = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\enfants';

freqVicon = 100;
fc_filter = 6;

% Dossier de sortie
output_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Fig\Smoothness_Spectra';
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

% Filtre passe-bas
[b, a] = butter(2, fc_filter/(freqVicon/2), 'low');

% Couleurs par surface (cohérentes avec tes autres scripts)
color_map = containers.Map( ...
    {'Plat', 'Medium', 'High'}, ...
    {[0.2 0.4 0.8], [0.2 0.7 0.3], [0.8 0.2 0.2]});  % Bleu, Vert, Rouge

%% === PARAMÈTRES SPARC (IDENTIQUES À TON CODE) ===
ampThreshold = 0.05;
fMaxHz       = 6;
zeroPadIdx   = 4;

%% === STRUCTURE POUR STOCKER LES SPECTRES ===
spectral_data = struct();

fprintf('🔄 Traitement spectral pour %s - Essai %d\n\n', sujet_id, essai_cible);

%% === BOUCLE DE TRAITEMENT PAR SURFACE ===
for surf_idx = 1:length(surfaces)
    surface = surfaces{surf_idx};
    
    filename = sprintf('%s_%s_%02d.c3d', sujet_id, surface, essai_cible);
    c3d_path = fullfile(base_dir, filename);
    
    if ~isfile(c3d_path)
        fprintf('⚠️  Fichier manquant : %s\n', filename);
        continue;
    end
    
    fprintf('📂 Traitement : %s\n', filename);
    
    try
        % === LECTURE C3D ===
        data = btkReadAcquisition(c3d_path);
        markers = btkGetMarkers(data);
        
        % === COM PELVIEN ===
        COM = calculate_pelvic_COM(markers);
        
        % === FENÊTRE D'ANALYSE : du 1er au dernier HS ===
        n_frames = size(COM, 1);
        fsrt_frame = btkGetFirstFrame(data);
        events = btkGetEvents(data);
        
        % HS gauche
        if isfield(events, 'Left_Foot_Strike')
            Left_HS_frames = round(events.Left_Foot_Strike * freqVicon - fsrt_frame + 1);
        else
            Left_HS_frames = [];
        end
        
        % HS droite
        if isfield(events, 'Right_Foot_Strike')
            Right_HS_frames = round(events.Right_Foot_Strike * freqVicon - fsrt_frame + 1);
        else
            Right_HS_frames = [];
        end
        
        HS_all = sort([Left_HS_frames(:); Right_HS_frames(:)]);
        clipHS = @(v) max(1, min(v, n_frames));
        HS_all = clipHS(HS_all);
        
        if numel(HS_all) >= 2
            start_frame = HS_all(1);
            end_frame = HS_all(end);
        else
            start_frame = 1;
            end_frame = n_frames;
        end
        
        if end_frame <= start_frame
            start_frame = 1;
            end_frame = n_frames;
        end
        
        fprintf('   📏 Analyse entre frames %d et %d\n', start_frame, end_frame);
        
        % === FILTRAGE POSITION ===
        COM_filt = filtfilt(b, a, COM);
        
        % === VITESSES ===
        vel_COM = calculate_velocities(COM_filt, freqVicon);
        
        % === MAGNITUDE 3D DE LA VITESSE ===
        v_mag = sqrt(sum(vel_COM(start_frame:end_frame, :).^2, 2));
        speed = v_mag(:)';  % 1xN (format requis pour SPARC)
        
        % === CALCUL DU SPECTRE (MÉTHODE SPARC) ===
        N = length(speed);
        Nfft = 2^(ceil(log2(N)) + zeroPadIdx);
        
        % FFT
        V = abs(fft(speed, Nfft));
        V = V(1:Nfft/2+1);  % Garder seulement les fréquences positives
        
        % Vecteur de fréquences
        freqs = (0:Nfft/2)' * (freqVicon / Nfft);
        
        % Normalisation par DC component (comme dans SPARC)
        V0 = max(V(1), 1e-10);
        Vn = V / V0;
        
        % Calcul de SPARC pour référence
        omega = 2 * pi * freqs;
        omega_c = 2 * pi * fMaxHz;
        idx_sparc = find(omega <= omega_c);
        
        if length(idx_sparc) >= 2
            omega_seg = omega(idx_sparc);
            V_seg = Vn(idx_sparc);
            d_omega = diff(omega_seg);
            d_V = diff(V_seg);
            arc = sum(sqrt((d_omega / omega_c).^2 + d_V.^2));
            sparc_value = -arc;
        else
            sparc_value = NaN;
        end
        
        % === STOCKAGE ===
        spectral_data.(surface).freqs = freqs;
        spectral_data.(surface).V_raw = V;           % Spectre non-normalisé
        spectral_data.(surface).V_norm = Vn;         % Spectre normalisé
        spectral_data.(surface).speed = speed;        % Signal temporel
        spectral_data.(surface).duration_sec = (end_frame - start_frame + 1) / freqVicon;
        spectral_data.(surface).sparc = sparc_value;
        spectral_data.(surface).n_steps = numel(HS_all);
        
        fprintf('   ✅ SPARC = %.3f, Durée = %.1f sec, N_steps = %d\n', ...
                sparc_value, spectral_data.(surface).duration_sec, spectral_data.(surface).n_steps);
        
        btkCloseAcquisition(data);
        
    catch ME
        fprintf('❌ Erreur : %s\n', ME.message);
    end
end

%% ================================================================
%  FIGURE 1 : SPECTRES NON-NORMALISÉS (PSD) - STYLE ARTICLE
% ================================================================
fig1 = figure('Name', 'Spectres non-normalisés - COM Magnitude', ...
              'Position', [100, 100, 1400, 500], 'Color', 'w');

for surf_idx = 1:length(surfaces)
    surface = surfaces{surf_idx};
    
    if ~isfield(spectral_data, surface)
        continue;
    end
    
    subplot(1, 3, surf_idx);
    hold on; box on; grid on;
    
    freqs = spectral_data.(surface).freqs;
    V_raw = spectral_data.(surface).V_raw;
    
    % Tracer le spectre brut
    plot(freqs, V_raw, 'Color', color_map(surface), 'LineWidth', 2);
    
    % Marquer les pics principaux (cadence et harmoniques)
    [pks, locs] = findpeaks(V_raw, 'MinPeakHeight', max(V_raw)*0.1, 'NPeaks', 4, 'SortStr', 'descend');
    plot(freqs(locs), pks, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    % Annotations des fréquences de pic
    for i = 1:min(2, length(locs))
        text(freqs(locs(i)), pks(i)*1.1, sprintf('%.2f Hz', freqs(locs(i))), ...
             'FontSize', 9, 'HorizontalAlignment', 'center');
    end
    
    % Ligne verticale à 6 Hz (limite SPARC)
    xline(6, '--k', 'LineWidth', 1.5, 'Alpha', 0.7, 'Label', '6 Hz');
    
    xlabel('Fréquence (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Magnitude (ua)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('%s\nSPARC = %.3f, %d steps', surface, ...
                  spectral_data.(surface).sparc, spectral_data.(surface).n_steps), ...
          'FontSize', 12, 'FontWeight', 'bold');
    
    xlim([0 12]);
    ylim([0 inf]);
    
    hold off;
end

sgtitle(sprintf('Spectres Non-normalisés - %s - Essai %d - COM Magnitude', sujet_id, essai_cible), ...
        'FontSize', 14, 'FontWeight', 'bold');

saveas(fig1, fullfile(output_dir, sprintf('%s_Essai%02d_Spectra_Raw.png', sujet_id, essai_cible)));
fprintf('✅ Figure 1 sauvegardée\n');

%% ================================================================
%  FIGURE 2 : SPECTRES NORMALISÉS (STYLE BECK ET AL. 2018)
% ================================================================
fig2 = figure('Name', 'Spectres normalisés - COM Magnitude', ...
              'Position', [150, 150, 1400, 500], 'Color', 'w');

for surf_idx = 1:length(surfaces)
    surface = surfaces{surf_idx};
    
    if ~isfield(spectral_data, surface)
        continue;
    end
    
    subplot(1, 3, surf_idx);
    hold on; box on; grid on;
    
    freqs = spectral_data.(surface).freqs;
    V_norm = spectral_data.(surface).V_norm;
    
    % Tracer le spectre normalisé
    plot(freqs, V_norm, 'Color', color_map(surface), 'LineWidth', 2.5);
    
    % Zone SPARC (0-6 Hz) - remplissage sous la courbe
    idx_sparc = freqs <= 6;
    freqs_sparc = freqs(idx_sparc);
    V_sparc = V_norm(idx_sparc);
    
    % S'assurer que ce sont des vecteurs colonnes
    freqs_sparc = freqs_sparc(:);
    V_sparc = V_sparc(:);
    
    fill([freqs_sparc; flipud(freqs_sparc)], ...
         [V_sparc; zeros(size(V_sparc))], ...
         color_map(surface), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    
    % Ligne verticale à 6 Hz
    xline(6, '--k', 'LineWidth', 1.5, 'Alpha', 0.7, 'Label', '6 Hz cutoff');
    
    % Ligne horizontale au seuil d'amplitude (0.05)
    yline(0.05, ':k', 'LineWidth', 1, 'Alpha', 0.5, 'Label', 'Amp. threshold');
    
    xlabel('Fréquence (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Magnitude normalisée', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('%s\nSPARC = %.3f', surface, spectral_data.(surface).sparc), ...
          'FontSize', 12, 'FontWeight', 'bold');
    
    xlim([0 12]);
    ylim([0 1.2]);
    
    hold off;
end

sgtitle(sprintf('Spectres Normalisés - %s - Essai %d - COM Magnitude', sujet_id, essai_cible), ...
        'FontSize', 14, 'FontWeight', 'bold');

saveas(fig2, fullfile(output_dir, sprintf('%s_Essai%02d_Spectra_Normalized.png', sujet_id, essai_cible)));
fprintf('✅ Figure 2 sauvegardée\n');

%% ================================================================
%  FIGURE 3 : SUPERPOSITION DES 3 CONDITIONS (NON-NORMALISÉ)
% ================================================================
fig3 = figure('Name', 'Comparaison Spectres - 3 Surfaces', ...
              'Position', [200, 200, 800, 600], 'Color', 'w');
hold on; box on; grid on;

leg_entries = {};

for surf_idx = 1:length(surfaces)
    surface = surfaces{surf_idx};
    
    if ~isfield(spectral_data, surface)
        continue;
    end
    
    freqs = spectral_data.(surface).freqs;
    V_raw = spectral_data.(surface).V_raw;
    
    % Tracer avec transparence
    plot(freqs, V_raw, 'Color', color_map(surface), ...
         'LineWidth', 2.5, 'DisplayName', surface);
    
    leg_entries{end+1} = sprintf('%s (SPARC=%.3f)', surface, spectral_data.(surface).sparc);
end

% Ligne verticale à 6 Hz
xline(6, '--k', 'LineWidth', 2, 'Alpha', 0.7, 'Label', '6 Hz');

xlabel('Fréquence (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Magnitude (ua)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Comparaison Spectres - %s - Essai %d\nCOM Magnitude (Non-normalisé)', ...
              sujet_id, essai_cible), ...
      'FontSize', 13, 'FontWeight', 'bold');

xlim([0 12]);
ylim([0 inf]);

legend(leg_entries, 'Location', 'northeast', 'FontSize', 10);

hold off;

saveas(fig3, fullfile(output_dir, sprintf('%s_Essai%02d_Comparison_Raw.png', sujet_id, essai_cible)));
fprintf('✅ Figure 3 sauvegardée\n');

%% ================================================================
%  FIGURE 4 : SUPERPOSITION DES 3 CONDITIONS (NORMALISÉ)
% ================================================================
fig4 = figure('Name', 'Comparaison Spectres Normalisés - 3 Surfaces', ...
              'Position', [250, 250, 800, 600], 'Color', 'w');
hold on; box on; grid on;

leg_entries = {};

for surf_idx = 1:length(surfaces)
    surface = surfaces{surf_idx};
    
    if ~isfield(spectral_data, surface)
        continue;
    end
    
    freqs = spectral_data.(surface).freqs;
    V_norm = spectral_data.(surface).V_norm;
    
    % Tracer
    plot(freqs, V_norm, 'Color', color_map(surface), ...
         'LineWidth', 2.5, 'DisplayName', surface);
    
    leg_entries{end+1} = sprintf('%s (SPARC=%.3f)', surface, spectral_data.(surface).sparc);
end

% Ligne verticale à 6 Hz
xline(6, '--k', 'LineWidth', 2, 'Alpha', 0.7, 'Label', '6 Hz');

% Ligne horizontale au seuil
yline(0.05, ':k', 'LineWidth', 1.5, 'Alpha', 0.5, 'Label', 'Threshold 0.05');

xlabel('Fréquence (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Magnitude normalisée', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Comparaison Spectres Normalisés - %s - Essai %d\nCOM Magnitude', ...
              sujet_id, essai_cible), ...
      'FontSize', 13, 'FontWeight', 'bold');

xlim([0 12]);
ylim([0 1.2]);

legend(leg_entries, 'Location', 'northeast', 'FontSize', 10);

hold off;

saveas(fig4, fullfile(output_dir, sprintf('%s_Essai%02d_Comparison_Normalized.png', sujet_id, essai_cible)));
fprintf('✅ Figure 4 sauvegardée\n');

%% ================================================================
%  TABLEAU RÉCAPITULATIF
% ================================================================
summary_table = table();
row_idx = 1;

for surf_idx = 1:length(surfaces)
    surface = surfaces{surf_idx};
    
    if ~isfield(spectral_data, surface)
        continue;
    end
    
    summary_table.Surface{row_idx} = surface;
    summary_table.SPARC(row_idx) = spectral_data.(surface).sparc;
    summary_table.Duration_sec(row_idx) = spectral_data.(surface).duration_sec;
    summary_table.N_Steps(row_idx) = spectral_data.(surface).n_steps;
    
    % Trouver la fréquence du pic principal
    [~, idx_max] = max(spectral_data.(surface).V_raw);
    summary_table.Peak_Freq_Hz(row_idx) = spectral_data.(surface).freqs(idx_max);
    
    row_idx = row_idx + 1;
end

disp('=== RÉSUMÉ SPECTRAL ===');
disp(summary_table);

writetable(summary_table, fullfile(output_dir, sprintf('%s_Essai%02d_Summary.csv', sujet_id, essai_cible)));
fprintf('✅ Tableau récapitulatif sauvegardé\n');

fprintf('\n🎯 TERMINÉ! 4 figures générées dans : %s\n', output_dir);

%% ============================== FONCTIONS ==============================

function COM = calculate_pelvic_COM(markers)
    pelvis_markers = {};
    needed = {'LPSI','RPSI','LASI','RASI'};

    for i = 1:length(needed)
        if isfield(markers, needed{i})
            mk = markers.(needed{i});
            if ~isempty(mk) && size(mk,2) == 3
                pelvis_markers{end+1} = mk;
            end
        end
    end

    if length(pelvis_markers) < 2
        error('Not enough pelvic markers (minimum 2 requis)');
    end

    L = min(cellfun(@(x) size(x,1), pelvis_markers));
    for i = 1:length(pelvis_markers)
        pelvis_markers{i} = pelvis_markers{i}(1:L,:);
    end

    stack = cat(3, pelvis_markers{:});
    COM   = mean(stack, 3);
end

function vel = calculate_velocities(pos, fs)
    vel = gradient(pos, 1/fs);
end