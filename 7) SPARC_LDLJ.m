%% === CALCUL SPARC/LDLJ - COM vs STERNUM ===
% Segmentation de chaque essai du 1er au dernier HS
% Recommandation des calculs par Balasubramanian et al. 2015

clear; clc; close all;

%% === CHEMINS ET PARAMÈTRES ===
addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\Scripts\btk'));    
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));

sujet_id  = 'CTL_01';
surfaces  = {'Plat', 'Medium', 'High'};
essais    = 1:10; % Attention entre les groupes !
base_dir  = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\enfants'; % Là où je vais chercher mes data

freqVicon = 100; % La frequence d'acquisition du Vicon
fc_filter = 6;

output_csv = sprintf('C:\\Users\\silve\\Desktop\\DOCTORAT\\UNIV MONTREAL\\TRAVAUX-THESE\\Surfaces_Irregulieres\\Datas\\Script\\gaitAnalysisGUI\\result\\Smoothness\\Smoothness_TrialBased_%s.csv', sujet_id);
output_mat = sprintf('C:\\Users\\silve\\Desktop\\DOCTORAT\\UNIV MONTREAL\\TRAVAUX-THESE\\Surfaces_Irregulieres\\Datas\\Script\\gaitAnalysisGUI\\result\\Smoothness\\Smoothness_TrialBased_%s.mat', sujet_id);

%% === INITIALISATION ===
results    = table();
log_errors = {};
fprintf('🔄 Analyse de fluidité PAR ESSAI COMPLET pour %s...\n\n', sujet_id);

% Filtre passe-bas
[b, a] = butter(2, fc_filter/(freqVicon/2), 'low');

%% === FLAG POUR VISUALISATION ===
PLOT_SPECTRES = true;  % Mettre à false pour désactiver les plots
MAX_PLOTS = 3;         % Nombre maximum de spectres à afficher
plot_counter = 0;

% === BOUCLE DE TRAITEMENT ===
for surf_idx = 1:length(surfaces)
    surface = surfaces{surf_idx};
    
    for essai = essais
        filename = sprintf('%s_%s_%02d.c3d', sujet_id, surface, essai);
        c3d_path = fullfile(base_dir, filename);

        if ~isfile(c3d_path)
            msg = sprintf('❌ Fichier manquant : %s', filename);
            log_errors{end+1} = msg;
            fprintf('%s\n', msg);
            continue;
        end
        
        fprintf('📂 Traitement : %s\n', filename);
        
        try
            % Lecture C3D
            data    = btkReadAcquisition(c3d_path);
            markers = btkGetMarkers(data);

            % === COM PELVIEN ===
            try
                COM = calculate_pelvic_COM(markers);   % [N x 3]
            catch
                msg = sprintf('❌ %s : COM impossible (markers manquants ou invalides)', filename);
                log_errors{end+1} = msg;
                fprintf('%s\n', msg);
                btkCloseAcquisition(data);
                continue;
            end

            % === STERNUM ===
            if isfield(markers,'STRN')
                STERN = markers.STRN;     % [N x 3] Plug-in Gait sternum
            else
                msg = sprintf('⚠️ %s : marqueur STRN absent, indices STERN mis à NaN', filename);
                log_errors{end+1} = msg;
                fprintf('%s\n', msg);
                STERN = NaN(size(COM));   % même taille, rempli de NaN
            end

            % === FENÊTRE D'ANALYSE : du 1er HS au dernier HS (toutes jambes confondues) ===
            n_frames = size(COM, 1);
            fsrt_frame = btkGetFirstFrame(data);
            events     = btkGetEvents(data);

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

            % Clip de sécurité dans [1, n_frames]
            clipHS = @(v) max(1, min(v, n_frames));
            HS_all = clipHS(HS_all);

            if numel(HS_all) >= 2
                start_frame = HS_all(1);
                end_frame   = HS_all(end);
            else
                % Fallback : si pas (ou pas assez) d'événements HS -> essai complet
                start_frame = 1;
                end_frame   = n_frames;
            end

            % Sécurité finale : éviter fenêtre inversée / trop courte
            if end_frame <= start_frame
                start_frame = 1;
                end_frame   = n_frames;
            end

            n_cycles = NaN;

            fprintf('   📏 Analyse entre 1er HS et dernier HS (frames %d à %d)\n', start_frame, end_frame);

            % === FILTRAGE POSITION ===
            COM_filt   = filtfilt(b, a, COM);
            STERN_filt = filtfilt(b, a, STERN);

            % === VITESSES ===
            vel_COM   = calculate_velocities(COM_filt,   freqVicon);  
            vel_STERN = calculate_velocities(STERN_filt, freqVicon);  % (NaN si STRN absent)

            % === INDICES DE FLUIDITÉ ===
            smooth_COM   = calculate_smoothness_indices(vel_COM,   start_frame, end_frame, freqVicon);
            smooth_STERN = calculate_smoothness_indices(vel_STERN, start_frame, end_frame, freqVicon);

            % Vérification scalaires
            fields_COM = fieldnames(smooth_COM);
            for f = 1:length(fields_COM)
                if ~isscalar(smooth_COM.(fields_COM{f}))
                    msg = sprintf('⚠️ %s : COM Champ %s non scalaire → NaN', filename, fields_COM{f});
                    log_errors{end+1} = msg;
                    fprintf('%s\n', msg);
                    smooth_COM.(fields_COM{f}) = NaN;
                end
            end

            fields_STERN = fieldnames(smooth_STERN);
            for f = 1:length(fields_STERN)
                if ~isscalar(smooth_STERN.(fields_STERN{f}))
                    msg = sprintf('⚠️ %s : STERN Champ %s non scalaire → NaN', filename, fields_STERN{f});
                    log_errors{end+1} = msg;
                    fprintf('%s\n', msg);
                    smooth_STERN.(fields_STERN{f}) = NaN;
                end
            end

            % % === VISUALISATION POUR LES PREMIERS ESSAIS ===
            % if PLOT_SPECTRES && plot_counter < MAX_PLOTS
            %     % Extraire la vitesse pour une direction (ex: AP)
            %     v_seg = vel_COM(start_frame:end_frame, 1); % 1 = AP
            % 
            %     % Calculer SPARC avec visualisation
            %     titleStr = sprintf('%s - %s - Essai %02d - COM AP', ...
            %                        sujet_id, surface, essai);
            %     sparc_test = compute_SPARC(abs(v_seg), freqVicon, 'plot', true, 'title', titleStr);
            %     plot_counter = plot_counter + 1;
            % end

            % === STOCKAGE ===
            new_row = struct();
            new_row.Sujet            = {sujet_id};
            new_row.Surface          = {surface};
            new_row.Essai            = essai;
            new_row.N_Cycles         = n_cycles;
            new_row.Start_Frame      = start_frame;
            new_row.End_Frame        = end_frame;
            new_row.Duration_frames  = end_frame - start_frame + 1;
            new_row.Duration_sec     = new_row.Duration_frames / freqVicon;

            % Ajout COM_*
            for f = 1:length(fields_COM)
                fname = fields_COM{f};
                new_row.(sprintf('COM_%s', fname)) = smooth_COM.(fname);
            end

            % Ajout STERN_*
            for f = 1:length(fields_STERN)
                fname = fields_STERN{f};
                new_row.(sprintf('STERN_%s', fname)) = smooth_STERN.(fname);
            end

            results = [results; struct2table(new_row)];

            % Fermeture BTK
            btkCloseAcquisition(data);
            
            fprintf('   ✅ Essai traité : durée %.1f sec\n\n', new_row.Duration_sec);
            
        catch ME
            msg = sprintf('❌ %s : %s', filename, ME.message);
            log_errors{end+1} = msg;
            fprintf('%s\n\n', msg);
        end
    end
end

%% === SAUVEGARDE ===
fprintf('💾 Sauvegarde...\n');

% .CSV
writetable(results, output_csv);
% .MAT
save(output_mat, 'results', 'log_errors');

nTrials = height(results);
fprintf('📊 Total essais traités : %d\n', nTrials);

%% === FIGURES COMPARATIVES COM vs STERN ===

% Familles d'indices et directions
families   = {'SPARC', 'LDLJ'};
directions = {'Magnitude', 'AP', 'ML', 'V'};

metrics_to_plot = {};

for f = 1:numel(families)
    for d = 1:numel(directions)
        metrics_to_plot{end+1} = sprintf('%s_%s', families{f}, directions{d});
    end
end

% Boucle sur toutes les métriques disponibles
for m = 1:numel(metrics_to_plot)
    metricName = metrics_to_plot{m};
    % On trace seulement si les colonnes existent vraiment
    comField   = ['COM_'   metricName];
    sternField = ['STERN_' metricName];

    if ismember(comField, results.Properties.VariableNames) && ...
       ismember(sternField, results.Properties.VariableNames)
        plot_COM_vs_STERN(results, metricName);
    else
        fprintf('⏭️  Skip %s (champs %s ou %s absents)\n', ...
            metricName, comField, sternField);
    end
end

%% === COMPARAISON DES SPECTRES PAR SURFACE ===
fprintf('\n\n🎨 === ANALYSE SPECTRALE COMPARATIVE ===\n');

% Appel de la fonction de comparaison
compare_spectres_surfaces(base_dir, sujet_id, surfaces, essais, freqVicon, fc_filter);

fprintf('📊 Analyse spectrale terminée!\n\n');

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

function smoothness = calculate_smoothness_indices(vel, i1, i2, fs)
    % Si vel contient des NaN (ex STRN absent) => tout à NaN
    if all(isnan(vel(:)))
        smoothness.SPARC_AP        = NaN;
        smoothness.SPARC_ML        = NaN;
        smoothness.SPARC_V         = NaN;
        smoothness.SPARC_Magnitude = NaN;
        smoothness.LDLJ_AP         = NaN;
        smoothness.LDLJ_ML         = NaN;
        smoothness.LDLJ_V          = NaN;
        smoothness.LDLJ_Magnitude  = NaN;
        return;
    end

    % --- Paramètres SPARC ---
    Ts = 1/fs;
    ampThreshold = 0.05;
    fMaxHz       = 6;   % cohérent avec passe-bas à 6 Hz et analyse de marche
    zeroPadIdx   = 4;
    params       = [ampThreshold, fMaxHz, zeroPadIdx];

    directions = {'AP','ML','V'};

    % === SPARC et LDLJ par direction ===
    for d = 1:3
        v_seg = vel(i1:i2, d);

        % IMPORTANT: SpectralArcLength attend un vecteur vitesse (speed) en LIGNE 1xN
        % Tu utilises abs(...) pour avoir un speed positif (comme dans ton code actuel)
        speed = abs(v_seg(:))';  % 1xN

        smoothness.(sprintf('SPARC_%s', directions{d})) = SpectralArcLength(speed, Ts, params);
        smoothness.(sprintf('LDLJ_%s',  directions{d})) = compute_LDLJ(v_seg, fs);
    end

    % === Magnitude 3D ===
    v_mag = sqrt(sum(vel(i1:i2,:).^2, 2));
    speed3D = v_mag(:)';  % 1xN

    smoothness.SPARC_Magnitude = SpectralArcLength(speed3D, Ts, params);
    smoothness.LDLJ_Magnitude  = compute_LDLJ(v_mag, fs);
end

% function sparc = compute_SPARC(velocity, fs, varargin)
% % compute_SPARC - Calcule le SPARC avec option de visualisation du spectre
% %
% % Syntaxe:
% %   sparc = compute_SPARC(velocity, fs)
% %   sparc = compute_SPARC(velocity, fs, 'plot', true)
% %   sparc = compute_SPARC(velocity, fs, 'plot', true, 'title', 'Mon titre')
% %
% % Paramètres:
% %   velocity : vecteur de vitesse
% %   fs       : fréquence d'échantillonnage (Hz)
% %   'plot'   : true/false pour afficher le spectre (défaut: false)
% %   'title'  : titre personnalisé pour la figure
% 
%     % Parsing des arguments optionnels
%     p = inputParser;
%     addParameter(p, 'plot', false, @islogical);
%     addParameter(p, 'title', '', @ischar);
%     parse(p, varargin{:});
% 
%     doPlot = p.Results.plot;
%     plotTitle = p.Results.title;
% 
%     % Vérifications initiales
%     if length(velocity) < 20 || all(isnan(velocity))
%         sparc = NaN;
%         return;
%     end
% 
%     velocity = velocity(:);
% 
%     % Etapes basées sur Balasubramanian et al. 2012
% 
%     % ÉTAPE 1: FFT avec zero-padding (padlevel = 4)
%     N = length(velocity);
%     K = 2^(nextpow2(N) + 4);  
%     V = abs(fft(velocity, K));
%     V = V(1:K/2+1);
% 
%     % ÉTAPE 2: Normalisation par DC
%     V0 = max(V(1), 1e-10);
%     Vn = V / V0; 
%     %Vn = V / max(V); % test normalisation par max du spectre (pas recommandé)
% 
%     % ÉTAPE 3: Arc length sur 0-6 Hz
%     freqs = (0:K/2)' * (fs / K);
%     f_c = 6;
%     omega = 2 * pi * freqs;
%     omega_c = 2 * pi * f_c;  % 6 Hz
%     idx = find(omega <= omega_c);
% 
%     if length(idx) < 2
%         sparc = NaN;
%         return;
%     end
% 
%     omega_seg = omega(idx);
%     V_seg = Vn(idx);
%     d_omega = diff(omega_seg);
%     d_V = diff(V_seg);
%     arc = sum(sqrt((d_omega / omega_c).^2 + d_V.^2));
%     sparc = -arc;
% 
%     % VISUALISATION
%     if doPlot
%         figure('Name', 'Analyse spectrale SPARC', 'Color', 'w');
% 
%         % Subplot 1: Signal temporel
%         subplot(3,1,1);
%         t = (0:length(velocity)-1) / fs;
%         plot(t, velocity, 'b-', 'LineWidth', 1.5);
%         grid on; box on;
%         xlabel('Temps (s)');
%         ylabel('Vitesse (m/s)');
%         title('Signal temporel');
% 
%         % Subplot 2: Spectre complet (0 à Nyquist)
%         subplot(3,1,2);
%         plot(freqs, Vn, 'k-', 'LineWidth', 1);
%         hold on;
%         xline(6, 'r--', 'LineWidth', 2, 'Label', '6 Hz (cutoff)');
%         hold off;
%         grid on; box on;
%         xlabel('Fréquence (Hz)');
%         ylabel('Magnitude normalisée');
%         title('Spectre fréquentiel complet');
%         xlim([0 min(50, fs/2)]);
% 
%         % Subplot 3: Spectre dans la bande SPARC (0-6 Hz) avec arc-length
%         subplot(3,1,3);
%         freq_seg = omega_seg / (2*pi);
%         plot(freq_seg, V_seg, 'b-', 'LineWidth', 2);
%         hold on;
% 
%         % Visualisation de l'arc-length
%         plot(freq_seg, V_seg, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
% 
%         hold off;
%         grid on; box on;
%         xlabel('Fréquence (Hz)');
%         ylabel('Magnitude normalisée');
%         title(sprintf('Bande SPARC (0-6 Hz) | SPARC = %.3f | Arc-length = %.3f', sparc, arc));
%         xlim([0 6]);
% 
%         % Titre global
%         if ~isempty(plotTitle)
%             sgtitle(plotTitle, 'FontSize', 12, 'FontWeight', 'bold');
%         end
%     end
% end

function ldlj = compute_LDLJ(velocity, fs) 
    if length(velocity) < 10 || all(isnan(velocity))
        ldlj = NaN;
        return;
    end
    
    velocity = velocity(:);
    dt       = 1 / fs;
    
    acc  = gradient(velocity, dt);
    jerk = gradient(acc, dt);
    
    T      = (length(velocity) - 1) * dt;
    v_peak = max(abs(velocity));
    
    if v_peak < 1e-8 || T < 1e-8
        ldlj = NaN;
        return;
    end
    
    integral_j = sum(jerk.^2) * dt;
    D          = (T^5 / v_peak^2) * integral_j;
    
    ldlj = -log(D + eps);
end

% LOCAL FUNCTION FROM GITHUB: https://github.com/siva82kb/smoothness/blob/master/matlab/SpectralArcLength.m
function S = SpectralArcLength( speed, Ts, parameters )
% SPECTRALARCLENGTH computes the smoothness of the give movement speed
% profile using the spectral arc length method.
% The this function takes three inputs and provides one output.
% Inputs: { speed, Ts*, parameters* } ('*' optional parameters)
%         SPEED: Speed it the speed profile of the movement. This is a 1xN
%         row vecotr. N is the total number of points in the speed profile.
%         The function assumes that the movement speed profile is already
%         filtered and segemented.
%
%         TS*: Sampling time in seconds. (DEFAULT VALUE = 0.01sec) NOTE: IF
%         YOUR DATA WAS NOT SAMPLED AT 100HZ, YOU MUST ENTER THE
%         APPROPRIATE SAMPLING TIME FOR ACCURATE RESULTS.
%
%         PARAMETERS*: This contains the parameters to be used spectral arc
%         lenght computation. This is a 1x2 column vector. This input
%         argument is option
%           - PARAMETER(1): The amplitude threshold to be used to choose
%           the cut-off frequency. The default value is chosen to be 0.05.
%           - PARAMETER(2): Maximum cut-off frequency for the spectral arc
%           length calcualtion. (DEFAULT VALUE = 10HZ) NOTE: 20Hz IS USED
%           TO REPRESENT THE MAXIMUM FREQUENCY COMPONENT OF A MOVEMENT.
%           THIS WAS CHOSEN TO COVER BOTH NORMAL AND ABNORMAL MOTOR
%           BEAHVIOUR. YOU CAN USE A VALUE LOWER THAN 20Hz IF YOU ARE AWARE
%           OF THE MAXIMUM FREQUENCY COMPONENT IN THE MOVEMENT OF INTEREST.
%           - PARAMETER(3): Zero padding index. This parameter controls the
%           resolution of the movement spectrum calculated from the speed
%           profile. (DEFAULT VALUE = 4). NOTE: IT IS NOT ADVISABLE TO USE
%           VALUES LOWER THAN 4.
%
% Outputs: { S }
%          S: This is smoothness of the given movement.
%
% For any queries about the method or the code, or if you come across any
% bugs in the code, feel free to contact me at siva82kb@gmail.com
% Sivakumar Balasubramanian. July 02, 2014.

% Check input arguments.
if nargin == 0
    disp('Error! Input at least the movement speed profile for the smoothness calculation.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
elseif nargin == 1
    % Default sampling time.
    Ts = 1/100; % 10ms.
elseif nargin == 2
    % Default parameters are use for the spectral arc length caclulations.
    parameters = [0.05, 10, 4];
end

% Check if the input argument are of the appropriate dimensions.
% Speed profile.
sz = size(speed);

% BUG FIX (minimum change):
% The original condition was inverted and rejected valid 1xN row vectors.
% This corrected condition throws an error ONLY when speed is NOT a 1xN row vector.
if ~((sz(1) == 1) && (sz(2) > 1))
    disp('Error! speed must be a row vector.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end

% Sampling time.
if ~isscalar(Ts)
    disp('Error! Ts must be a scalar.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end

% Parameters.
if length(parameters) ~= 3
    disp('Error! parameter is a vector with two elements.');
    fprintf('\n');
    help('SpectralArcLength');
    return;
end

% Calculate the spectrum of the speed profile.
N = length(speed);
Nfft = 2^(ceil(log2(N))+parameters(3));
speedSpectrum = abs(fft( speed, Nfft ));

% Normalize spectrum with respect to the DC component.
% NOTE: keep their original normalization (by max) to minimize changes.
% Also build freq with a guaranteed length matching Nfft to avoid index mismatch.
freq = (0:Nfft-1) * ((1/Ts)/Nfft);
speedSpectrum = speedSpectrum'/max(speedSpectrum);

% Choose the spectrum that is always above the amplitude threshold, and
% within the cut-off frequency.
inxFc = find((freq(1:end) <= parameters(2)) & ...
             (speedSpectrum(1:end) >= parameters(1)), 1, 'last');

% SAFETY: if nothing found (or numerical mismatch), clamp index
if isempty(inxFc)
    inxFc = length(speedSpectrum);
elseif inxFc > length(speedSpectrum)
    inxFc = length(speedSpectrum);
end

% Calculate the spectral arc length.
% 1. select the spectrum of interest.
speedSpectrum = speedSpectrum(1:inxFc);

% 2. Calculate the incremental arc lengths.
if inxFc < 2
    S = NaN;
    return;
end
dArcLengths = sqrt((1/(inxFc-1))^2 + (diff(speedSpectrum)).^2);

% 4. Compute movement smoothness.
S = -sum(dArcLengths);

return;
end

function plot_COM_vs_STERN(results, metricFieldBase)
% Crée une figure comparant COM vs STERN pour une métrique donnée
% metricFieldBase = 'SPARC_Magnitude', 'SPARC_AP', 'LDLJ_ML', etc.
    
    % Noms complets des champs
    comField   = ['COM_'   metricFieldBase];
    sternField = ['STERN_' metricFieldBase];

    % Vérification de la présence des champs
    if ~ismember(comField, results.Properties.VariableNames) || ...
       ~ismember(sternField, results.Properties.VariableNames)
        warning('Champs %s ou %s introuvables dans results.', comField, sternField);
        return;
    end

    if ismember('RowType', results.Properties.VariableNames)
        trial_idx = results.RowType == "Trial";
    else
        trial_idx = true(height(results),1);
    end

    res_trial = results(trial_idx, :);

    % Surfaces dans l'ordre d'apparition parmi les essais
    surfaces = unique(res_trial.Surface, 'stable');
    nSurf    = numel(surfaces);

    % Couleurs COM vs STERN
    colCOM   = [0.2 0.4 0.8];
    colSTERN = [0.9 0.5 0.1];

    mCOM    = nan(nSurf,1);
    mSTERN  = nan(nSurf,1);
    sdCOM   = nan(nSurf,1);
    sdSTERN = nan(nSurf,1);

    dataCOM   = cell(nSurf,1);
    dataSTERN = cell(nSurf,1);

    % Récupère les données par surface
    for s = 1:nSurf
        idx = strcmp(res_trial.Surface, surfaces{s});

        yCOM   = res_trial.(comField)(idx);
        ySTERN = res_trial.(sternField)(idx);

        dataCOM{s}   = yCOM;
        dataSTERN{s} = ySTERN;

        mCOM(s)    = mean(yCOM,   'omitnan');
        mSTERN(s)  = mean(ySTERN, 'omitnan');
        sdCOM(s)   = std(yCOM,   'omitnan');
        sdSTERN(s) = std(ySTERN, 'omitnan');
    end

    % Création de la figure
    figure('Name', ['COM vs STERN - ' metricFieldBase], 'Color', 'w');
    hold on; box on;

    x        = 1:nSurf;
    barWidth = 0.35;

    % Barres des moyennes
    b1 = bar(x - barWidth/2, mCOM,   barWidth, 'FaceColor', colCOM,   'FaceAlpha', 0.6, 'EdgeColor', 'none');
    b2 = bar(x + barWidth/2, mSTERN, barWidth, 'FaceColor', colSTERN, 'FaceAlpha', 0.6, 'EdgeColor', 'none');

    % Barres d'erreur (±1 SD)
    for s = 1:nSurf
        if ~isnan(mCOM(s)) && ~isnan(sdCOM(s))
            plot([x(s)-barWidth/2 x(s)-barWidth/2], [mCOM(s)-sdCOM(s) mCOM(s)+sdCOM(s)], 'k-', 'LineWidth', 1);
        end
        if ~isnan(mSTERN(s)) && ~isnan(sdSTERN(s))
            plot([x(s)+barWidth/2 x(s)+barWidth/2], [mSTERN(s)-sdSTERN(s) mSTERN(s)+sdSTERN(s)], 'k-', 'LineWidth', 1);
        end
    end

    % Points individuels (essais)
    jitter = barWidth/3;

    for s = 1:nSurf
        yC = dataCOM{s};
        yS = dataSTERN{s};

        % On enlève les NaN pour éviter des points invisibles
        yC = yC(~isnan(yC));
        yS = yS(~isnan(yS));

        % COM – points
        if ~isempty(yC)
            xC = (x(s)-barWidth/2) + (rand(size(yC))-0.5) * jitter;
            scatter(xC, yC, 30, colCOM, 'filled', ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor','none');
        end

        % STERN – points
        if ~isempty(yS)
            xS = (x(s)+barWidth/2) + (rand(size(yS))-0.5) * jitter;
            scatter(xS, yS, 30, colSTERN, 'filled', ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor','none');
        end
    end

    % Axes & légendes
    set(gca, 'XTick', x, 'XTickLabel', surfaces, 'FontSize', 12);
    xlabel('Surface');
    ylabel(metricFieldBase, 'Interpreter','none');
    legend([b1 b2], {'COM', 'STERNUM'}, 'Location', 'best');
    title(sprintf('COM vs STERN - %s', metricFieldBase), 'Interpreter','none');

    hold off;
end

function compare_spectres_surfaces(base_dir, sujet_id, surfaces, essais, freqVicon, fc_filter)
% Compare les spectres fréquentiels entre les 3 surfaces
% Pour COM et STERN, directions AP, ML, V et Magnitude
%
% CORRECTION: Utilise une grille de fréquences commune pour éviter les
% erreurs de dimensions lors de la concaténation des spectres

    % Filtre passe-bas
    [b, a] = butter(2, fc_filter/(freqVicon/2), 'low');
    
    % Couleurs pour chaque surface
    colors = struct();
    colors.Plat = [0.2 0.6 0.2];    % Vert
    colors.Medium = [0.9 0.6 0.1];   % Orange
    colors.High = [0.8 0.2 0.2];     % Rouge
    
    % Stockage des spectres moyens par surface
    spectres_COM = struct();
    spectres_STERN = struct();
    directions = {'AP', 'ML', 'V', 'Magnitude'};
    
    % CORRECTION: Grille de fréquences commune pour interpolation (0 à 50 Hz par pas de 0.1 Hz)
    % Cela garantit que tous les spectres ont exactement la même longueur
    freqs_common = (0:0.1:50)';
    
    % Initialisation
    for surf_idx = 1:length(surfaces)
        surf = surfaces{surf_idx};
        for d = 1:length(directions)
            dir = directions{d};
            spectres_COM.(surf).(dir) = [];
            spectres_STERN.(surf).(dir) = [];
        end
    end
    
    % === COLLECTE DES DONNÉES ===
    fprintf('\n🔍 Collecte des spectres pour comparaison...\n');
    
    for surf_idx = 1:length(surfaces)
        surface = surfaces{surf_idx};
        fprintf('  📊 Surface: %s\n', surface);
        
        for essai = essais
            filename = sprintf('%s_%s_%02d.c3d', sujet_id, surface, essai);
            c3d_path = fullfile(base_dir, filename);
            
            if ~isfile(c3d_path)
                continue;
            end
            
            try
                % Lecture C3D
                data = btkReadAcquisition(c3d_path);
                markers = btkGetMarkers(data);
                
                % COM
                COM = calculate_pelvic_COM(markers);
                COM_filt = filtfilt(b, a, COM);
                vel_COM = calculate_velocities(COM_filt, freqVicon);
                
                % STERNUM
                if isfield(markers,'STRN')
                    STERN = markers.STRN;
                    STERN_filt = filtfilt(b, a, STERN);
                    vel_STERN = calculate_velocities(STERN_filt, freqVicon);
                else
                    vel_STERN = NaN(size(vel_COM));
                end
                
                % Calcul des spectres pour chaque direction
                for d = 1:3
                    dir = directions{d};
                    
                    % COM
                    [freqs, Vn_COM] = compute_spectrum(abs(vel_COM(:,d)), freqVicon);
                    % CORRECTION: Interpolation sur la grille commune
                    Vn_COM_interp = interp1(freqs, Vn_COM, freqs_common, 'linear', 0);
                    spectres_COM.(surface).(dir) = [spectres_COM.(surface).(dir); Vn_COM_interp'];
                    
                    % STERN
                    if ~all(isnan(vel_STERN(:)))
                        [~, Vn_STERN] = compute_spectrum(abs(vel_STERN(:,d)), freqVicon);
                        % CORRECTION: Interpolation sur la grille commune
                        Vn_STERN_interp = interp1(freqs, Vn_STERN, freqs_common, 'linear', 0);
                        spectres_STERN.(surface).(dir) = [spectres_STERN.(surface).(dir); Vn_STERN_interp'];
                    end
                end
                
                % Magnitude
                v_mag_COM = sqrt(sum(vel_COM.^2, 2));
                [freqs, Vn_mag_COM] = compute_spectrum(v_mag_COM, freqVicon);
                % CORRECTION: Interpolation sur la grille commune
                Vn_mag_COM_interp = interp1(freqs, Vn_mag_COM, freqs_common, 'linear', 0);
                spectres_COM.(surface).Magnitude = [spectres_COM.(surface).Magnitude; Vn_mag_COM_interp'];
                
                if ~all(isnan(vel_STERN(:)))
                    v_mag_STERN = sqrt(sum(vel_STERN.^2, 2));
                    [~, Vn_mag_STERN] = compute_spectrum(v_mag_STERN, freqVicon);
                    % CORRECTION: Interpolation sur la grille commune
                    Vn_mag_STERN_interp = interp1(freqs, Vn_mag_STERN, freqs_common, 'linear', 0);
                    spectres_STERN.(surface).Magnitude = [spectres_STERN.(surface).Magnitude; Vn_mag_STERN_interp'];
                end
                
                btkCloseAcquisition(data);
                
            catch ME
                warning('Erreur lors du traitement de %s: %s', filename, ME.message);
                continue;
            end
        end
    end
    
    % === CRÉATION DES FIGURES ===
    fprintf('\n📈 Création des figures comparatives...\n');
    
    % Figure 1: COM - toutes directions
    plot_comparison_figure(spectres_COM, freqs_common, surfaces, colors, 'COM', directions);
    
    % Figure 2: STERN - toutes directions
    plot_comparison_figure(spectres_STERN, freqs_common, surfaces, colors, 'STERNUM', directions);
    
    % Figure 3: COM vs STERN par direction
    plot_COM_vs_STERN_spectres(spectres_COM, spectres_STERN, freqs_common, surfaces, colors, directions);
    
    fprintf('✅ Figures créées!\n');
end

%  === FONCTION: CALCUL DU SPECTRE ===
function [freqs, Vn] = compute_spectrum(velocity, fs)
% Calcule le spectre de puissance normalisé d'un signal de vitesse
% Retourne les fréquences et les magnitudes normalisées

    velocity = velocity(:);
    
    % FFT avec zero-padding (identique à compute_SPARC)
    N = length(velocity);
    K = 2^(nextpow2(N) + 4);
    V = abs(fft(velocity, K));
    V = V(1:K/2+1);
    
    % Normalisation par la composante DC
    V0 = max(V(1), 1e-10);
    Vn = V / V0;
    
    % Vecteur de fréquences
    freqs = (0:K/2)' * (fs / K);
end

% === FONCTION: FIGURE DE COMPARAISON ===
function plot_comparison_figure(spectres, freqs, surfaces, colors, label, directions)
% Crée une figure 2x2 comparant les spectres des 3 surfaces
% pour chaque direction (AP, ML, V, Magnitude)
    
    figure('Name', sprintf('Comparaison Spectres - %s', label), ...
           'Position', [100 100 1600 900], 'Color', 'w');
    
    for d = 1:length(directions)
        dir = directions{d};
        subplot(2, 2, d);
        hold on; box on; grid on;
        
        leg_entries = {};
        
        for surf_idx = 1:length(surfaces)
            surf = surfaces{surf_idx};
            
            if isempty(spectres.(surf).(dir))
                continue;
            end
            
            % Calcul de la moyenne et écart-type sur les essais
            mean_spectrum = mean(spectres.(surf).(dir), 1, 'omitnan');
            std_spectrum = std(spectres.(surf).(dir), 0, 1, 'omitnan');
            
            % Limiter à 0-12 Hz pour visualisation
            idx_plot = freqs <= 12;
            f_plot = freqs(idx_plot);
            m_plot = mean_spectrum(idx_plot);
            s_plot = std_spectrum(idx_plot);
            
            % Plot moyenne avec zone d'écart-type
            col = colors.(surf);
            
            % Zone d'écart-type (transparente)
            fill([f_plot; flipud(f_plot)], ...
                 [m_plot' + s_plot'; flipud(m_plot' - s_plot')], ...
                 col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            
            % Ligne moyenne (épaisse)
            plot(f_plot, m_plot, '-', 'Color', col, 'LineWidth', 2.5);
            
            leg_entries{end+1} = surf;
        end
        
        % Ligne verticale à 6 Hz (limite SPARC)
        xline(6, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
        
        xlabel('Fréquence (Hz)', 'FontSize', 11);
        ylabel('Magnitude normalisée', 'FontSize', 11);
        title(sprintf('%s - %s', label, dir), 'FontSize', 12, 'FontWeight', 'bold');
        xlim([0 12]);
        ylim([0 inf]);
        
        if ~isempty(leg_entries)
            legend(leg_entries, 'Location', 'northeast', 'FontSize', 10);
        end
        
        hold off;
    end
    
    sgtitle(sprintf('Comparaison des spectres fréquentiels - %s (moyenne ± SD)', label), ...
            'FontSize', 14, 'FontWeight', 'bold');
end

% === FONCTION: COM vs STERN PAR SURFACE ===
function plot_COM_vs_STERN_spectres(spectres_COM, spectres_STERN, freqs, surfaces, colors, directions)
% Crée une figure par surface comparant COM vs STERN
% pour les 4 directions (AP, ML, V, Magnitude)
    
    % Une figure par surface
    for surf_idx = 1:length(surfaces)
        surf = surfaces{surf_idx};
        
        figure('Name', sprintf('COM vs STERN - %s', surf), ...
               'Position', [100 100 1600 900], 'Color', 'w');
        
        col_COM = [0.2 0.4 0.8];
        col_STERN = [0.9 0.5 0.1];
        
        for d = 1:length(directions)
            dir = directions{d};
            subplot(2, 2, d);
            hold on; box on; grid on;
            
            % COM
            if ~isempty(spectres_COM.(surf).(dir))
                mean_COM = mean(spectres_COM.(surf).(dir), 1, 'omitnan');
                std_COM = std(spectres_COM.(surf).(dir), 0, 1, 'omitnan');
                
                idx_plot = freqs <= 12;
                f_plot = freqs(idx_plot);
                
                % Zone d'écart-type COM
                fill([f_plot; flipud(f_plot)], ...
                     [mean_COM(idx_plot)' + std_COM(idx_plot)'; ...
                      flipud(mean_COM(idx_plot)' - std_COM(idx_plot)')], ...
                     col_COM, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                
                plot(f_plot, mean_COM(idx_plot), '-', 'Color', col_COM, 'LineWidth', 2.5);
            end
            
            % STERN
            if ~isempty(spectres_STERN.(surf).(dir))
                mean_STERN = mean(spectres_STERN.(surf).(dir), 1, 'omitnan');
                std_STERN = std(spectres_STERN.(surf).(dir), 0, 1, 'omitnan');
                
                idx_plot = freqs <= 12;
                f_plot = freqs(idx_plot);
                
                % Zone d'écart-type STERN
                fill([f_plot; flipud(f_plot)], ...
                     [mean_STERN(idx_plot)' + std_STERN(idx_plot)'; ...
                      flipud(mean_STERN(idx_plot)' - std_STERN(idx_plot)')], ...
                     col_STERN, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                
                plot(f_plot, mean_STERN(idx_plot), '-', 'Color', col_STERN, 'LineWidth', 2.5);
            end
            
            % Ligne verticale à 20 Hz
            xline(6, 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
            
            xlabel('Fréquence (Hz)', 'FontSize', 11);
            ylabel('Magnitude normalisée', 'FontSize', 11);
            title(sprintf('%s', dir), 'FontSize', 12, 'FontWeight', 'bold');
            xlim([0 12]);
            ylim([0 inf]);
            legend({'COM', 'STERNUM'}, 'Location', 'northeast', 'FontSize', 10);
            
            hold off;
        end
        
        sgtitle(sprintf('COM vs STERNUM - Surface %s (moyenne ± SD)', surf), ...
                'FontSize', 14, 'FontWeight', 'bold');
    end
end