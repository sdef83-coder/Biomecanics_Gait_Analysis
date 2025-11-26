%% ================================================================
%  SCRIPT DE DIAGNOSTIC : COMPARAISON DES CALCULS DE MoS
%  Compare les valeurs de MoS_AP au Heel Strike entre :
%  - MOS.m (calcul direct)
%  - Visual_GlobalParticipant_MOS.m (avec interpolation)
%  - SpatioTemporal_Analysis.m (agrégation)
% ================================================================

clear; clc; close all;

%% === PARAMÈTRES À CONFIGURER ===
sujet_test = 'CTL_14';           % Participant à tester
surface_test = 'Plat';           % Surface à tester
essai_test = 1;                  % Numéro d'essai

% Chemins
base_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\jeunes_enfants';
mos_res_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\MoS';
addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\Scripts\btk'));
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL'))

% Charger L0
load('l0_participants.mat', 'l0_map');
if isKey(l0_map, sujet_test)
    L0_mm = l0_map(sujet_test) * 1000;
else
    error('Participant %s non trouvé dans l0_map', sujet_test);
end

%% === ÉTAPE 1 : CALCUL DIRECT (comme MOS.m) ===
fprintf('\n========== MÉTHODE 1 : CALCUL DIRECT (MOS.m) ==========\n');

filename = sprintf('%s_%s_%02d.c3d', sujet_test, surface_test, essai_test);
c3d_path = fullfile(base_dir, filename);

if ~isfile(c3d_path)
    error('Fichier C3D non trouvé : %s', c3d_path);
end

% Lecture du fichier
data = btkReadAcquisition(c3d_path);
markers = btkGetMarkers(data);
fsrt_frame = btkGetFirstFrame(data);
moments = btkGetEvents(data);
freqVicon = 100;

% Extraction des événements
if isfield(moments, 'Left_Foot_Strike')
    Left_HS_frames = round(moments.Left_Foot_Strike * freqVicon - fsrt_frame + 1);
else
    Left_HS_frames = [];
end

if isfield(moments, 'Left_Foot_Off')
    Left_TO_frames = round(moments.Left_Foot_Off * freqVicon - fsrt_frame + 1);
else
    Left_TO_frames = [];
end

% Sécurité
fnm = fieldnames(markers);
N = size(markers.(fnm{1}), 1);
clip = @(v) max(1, min(v, N-1));
Left_HS_frames = clip(Left_HS_frames);
Left_TO_frames = clip(Left_TO_frames);

fprintf('Cycles gauches détectés : %d\n', length(Left_HS_frames));

% Traitement du premier cycle
if ~isempty(Left_HS_frames)
    hs = Left_HS_frames(1);
    valid_TO = Left_TO_frames(Left_TO_frames > hs);
    
    if ~isempty(valid_TO)
        to = valid_TO(1);
        ho = estimate_heel_off_pair(hs, to, markers, 'L');
        
        fprintf('  Cycle 1 : HS=%d, HO=%d, TO=%d\n', hs, ho, to);
        
        % Calcul MoS avec filtrage 6 Hz (MOS.m)
        mos_method1 = calculate_MoS_method1(c3d_path, hs, ho, to, L0_mm);
        
        fprintf('\n--- RÉSULTATS MÉTHODE 1 (filtrage 6 Hz) ---\n');
        fprintf('  MoS_AP_HS (mm)  : %.3f\n', mos_method1.MoS_Heel_Strike_AP);
        fprintf('  MoS_AP_HS (%%L0) : %.3f\n', mos_method1.MoS_Heel_Strike_AP_P);
        fprintf('  MoS_ML_HS (mm)  : %.3f\n', mos_method1.MoS_Heel_Strike_ML);
        fprintf('  MoS_ML_HS (%%L0) : %.3f\n', mos_method1.MoS_Heel_Strike_ML_P);
    else
        error('Pas de Toe-Off valide pour ce cycle');
    end
else
    error('Aucun Heel Strike détecté');
end

%% === ÉTAPE 2 : AVEC INTERPOLATION (comme Visual_GlobalParticipant_MOS.m) ===
fprintf('\n========== MÉTHODE 2 : AVEC INTERPOLATION (Visual) ==========\n');

n_points = 100;
curve = extract_full_mos_curve_method2(base_dir, sujet_test, surface_test, essai_test, hs, ho, to, n_points, L0_mm);

fprintf('\n--- RÉSULTATS MÉTHODE 2 (interpolation) ---\n');
fprintf('  MoS_AP_HS (mm)  : %.3f (1er point interpolé)\n', curve.mos_ap(1));
fprintf('  MoS_AP_HS (%%L0) : %.3f\n', curve.mos_ap_P(1));
fprintf('  MoS_ML_HS (mm)  : %.3f\n', curve.mos_ml(1));
fprintf('  MoS_ML_HS (%%L0) : %.3f\n', curve.mos_ml_P(1));

%% === ÉTAPE 3 : VÉRIFICATION AVEC LE FICHIER .MAT SAUVEGARDÉ ===
fprintf('\n========== MÉTHODE 3 : FICHIER .MAT SAUVEGARDÉ ==========\n');

mat_file = fullfile(mos_res_dir, sprintf('MoS_results_%s.mat', sujet_test));

if isfile(mat_file)
    tmp = load(mat_file, 'MoS_data');
    results = tmp.MoS_data.results;
    
    % Filtrer par surface et premier cycle
    idx_surf = strcmp(results.Surface, surface_test) & results.Essai == essai_test & results.Cycle == 1;
    
    if any(idx_surf)
        row = results(idx_surf, :);
        
        fprintf('\n--- RÉSULTATS MÉTHODE 3 (.mat sauvegardé) ---\n');
        fprintf('  MoS_AP_HS (mm)  : %.3f\n', row.MoS_Heel_Strike_AP);
        fprintf('  MoS_AP_HS (%%L0) : %.3f\n', row.MoS_Heel_Strike_AP_P);
        fprintf('  MoS_ML_HS (mm)  : %.3f\n', row.MoS_Heel_Strike_ML);
        fprintf('  MoS_ML_HS (%%L0) : %.3f\n', row.MoS_Heel_Strike_ML_P);
        
        % Comparaison
        fprintf('\n========== COMPARAISON DES MÉTHODES ==========\n');
        fprintf('MoS_AP_HS (mm) :\n');
        fprintf('  Méthode 1 vs 3 : Δ = %.3f mm\n', abs(mos_method1.MoS_Heel_Strike_AP - row.MoS_Heel_Strike_AP));
        fprintf('  Méthode 2 vs 3 : Δ = %.3f mm\n', abs(curve.mos_ap(1) - row.MoS_Heel_Strike_AP));
        fprintf('  Méthode 1 vs 2 : Δ = %.3f mm\n', abs(mos_method1.MoS_Heel_Strike_AP - curve.mos_ap(1)));
        
        fprintf('\nMoS_ML_HS (mm) :\n');
        fprintf('  Méthode 1 vs 3 : Δ = %.3f mm\n', abs(mos_method1.MoS_Heel_Strike_ML - row.MoS_Heel_Strike_ML));
        fprintf('  Méthode 2 vs 3 : Δ = %.3f mm\n', abs(curve.mos_ml(1) - row.MoS_Heel_Strike_ML));
        fprintf('  Méthode 1 vs 2 : Δ = %.3f mm\n', abs(mos_method1.MoS_Heel_Strike_ML - curve.mos_ml(1)));
        
    else
        warning('Cycle non trouvé dans le fichier .mat');
    end
else
    warning('Fichier .mat non trouvé : %s', mat_file);
end

%% === ÉTAPE 4 : TEST AVEC DIFFÉRENTS FILTRAGES ===
fprintf('\n========== EFFET DU FILTRAGE ==========\n');

% Test avec filtrage 4 Hz (comme dans Visual)
mos_4hz = calculate_MoS_with_filter(c3d_path, hs, ho, to, L0_mm, 4, 4);
fprintf('Filtrage 4 Hz (ordre 4) :\n');
fprintf('  MoS_AP_HS : %.3f mm\n', mos_4hz.MoS_Heel_Strike_AP);

% Test avec filtrage 6 Hz ordre 2 (comme MOS.m)
mos_6hz = calculate_MoS_with_filter(c3d_path, hs, ho, to, L0_mm, 6, 2);
fprintf('Filtrage 6 Hz (ordre 2) :\n');
fprintf('  MoS_AP_HS : %.3f mm\n', mos_6hz.MoS_Heel_Strike_AP);

fprintf('\nDifférence due au filtrage : %.3f mm\n', abs(mos_4hz.MoS_Heel_Strike_AP - mos_6hz.MoS_Heel_Strike_AP));

%% === VISUALISATION ===
fprintf('\n========== GÉNÉRATION DES GRAPHIQUES ==========\n');

figure('Name', 'Diagnostic MoS', 'Position', [100 100 1400 800]);

% Subplot 1 : Comparaison des courbes AP
subplot(2,2,1); hold on; grid on;
x_norm = linspace(0, 100, n_points);
plot(x_norm, curve.mos_ap, 'b-', 'LineWidth', 2, 'DisplayName', 'Méthode 2 (interp)');
yline(mos_method1.MoS_Heel_Strike_AP, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Méthode 1 (direct) au HS');
if exist('row', 'var')
    yline(row.MoS_Heel_Strike_AP, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Méthode 3 (.mat)');
end
xlabel('% du cycle d''appui');
ylabel('MoS AP (mm)');
title('Comparaison MoS AP');
legend('Location', 'best');

% Subplot 2 : Comparaison des courbes ML
subplot(2,2,2); hold on; grid on;
plot(x_norm, curve.mos_ml, 'b-', 'LineWidth', 2, 'DisplayName', 'Méthode 2 (interp)');
yline(mos_method1.MoS_Heel_Strike_ML, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Méthode 1 (direct) au HS');
if exist('row', 'var')
    yline(row.MoS_Heel_Strike_ML, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Méthode 3 (.mat)');
end
xlabel('% du cycle d''appui');
ylabel('MoS ML (mm)');
title('Comparaison MoS ML');
legend('Location', 'best');

% Subplot 3 : Effet du filtrage sur COM_y
subplot(2,2,3); hold on; grid on;
data_raw = btkReadAcquisition(c3d_path);
markers_raw = btkGetMarkers(data_raw);
COM_y_raw = mean([markers_raw.LPSI(:,2), markers_raw.RPSI(:,2), markers_raw.LASI(:,2), markers_raw.RASI(:,2)], 2);

% Filtrage 6 Hz
[b6, a6] = butter(2, 6/50, 'low');
COM_y_6hz = filtfilt(b6, a6, COM_y_raw);

% Filtrage 4 Hz
[b4, a4] = butter(4, 4/50, 'low');
COM_y_4hz = filtfilt(b4, a4, COM_y_raw);

idx_plot = max(1, hs-50):min(length(COM_y_raw), to+50);
plot(idx_plot, COM_y_raw(idx_plot), 'k-', 'LineWidth', 0.8, 'DisplayName', 'Brut');
plot(idx_plot, COM_y_6hz(idx_plot), 'r-', 'LineWidth', 1.5, 'DisplayName', '6 Hz (ordre 2)');
plot(idx_plot, COM_y_4hz(idx_plot), 'b-', 'LineWidth', 1.5, 'DisplayName', '4 Hz (ordre 4)');
xline(hs, 'k--', 'HS'); xline(ho, 'm--', 'HO'); xline(to, 'g--', 'TO');
xlabel('Frame'); ylabel('COM_y (mm)');
title('Effet du filtrage sur COM_y');
legend('Location', 'best');

% Subplot 4 : Tableau récapitulatif
subplot(2,2,4);
axis off;

str_recap = {sprintf('PARTICIPANT : %s', sujet_test), ...
             sprintf('SURFACE : %s | ESSAI : %d', surface_test, essai_test), ...
             '', ...
             'MoS AP au HS (mm) :', ...
             sprintf('  Méthode 1 (direct 6Hz) : %.3f', mos_method1.MoS_Heel_Strike_AP), ...
             sprintf('  Méthode 2 (interp) : %.3f', curve.mos_ap(1)), ...
             ''};

if exist('row', 'var')
    str_recap{end+1} = sprintf('  Méthode 3 (.mat) : %.3f', row.MoS_Heel_Strike_AP);
    str_recap{end+1} = '';
    str_recap{end+1} = 'DIFFÉRENCES :';
    str_recap{end+1} = sprintf('  M1 vs M3 : %.3f mm', abs(mos_method1.MoS_Heel_Strike_AP - row.MoS_Heel_Strike_AP));
    str_recap{end+1} = sprintf('  M2 vs M3 : %.3f mm', abs(curve.mos_ap(1) - row.MoS_Heel_Strike_AP));
end

str_recap{end+1} = '';
str_recap{end+1} = sprintf('Effet filtrage 4Hz vs 6Hz : %.3f mm', abs(mos_4hz.MoS_Heel_Strike_AP - mos_6hz.MoS_Heel_Strike_AP));

text(0.1, 0.9, str_recap, 'VerticalAlignment', 'top', 'FontSize', 10, 'FontName', 'Courier');

fprintf('\n✅ DIAGNOSTIC TERMINÉ\n');

btkCloseAcquisition(data);
if exist('data_raw', 'var'), btkCloseAcquisition(data_raw); end

%% ========== FONCTIONS AUXILIAIRES ==========

function ho = estimate_heel_off_pair(hs, to, markers, side)
    heel_marker = [side 'HEE'];
    
    if ~isfield(markers, heel_marker) || ~(to > hs)
        ho = hs + 1;
        return
    end
    
    heel_z = markers.(heel_marker)(hs:to, 3);
    heel_z = fillmissing(heel_z, 'linear', 'EndValues', 'nearest');
    
    fs = 100;
    [bb, aa] = butter(2, 6/(fs/2), 'low');
    if numel(heel_z) > 6
        heel_z = filtfilt(bb, aa, heel_z);
    end
    
    heel_vel_z = diff(heel_z);
    stanceLen = to - hs;
    
    if stanceLen < 15
        ho = round(hs + 0.60*stanceLen);
        return
    end
    
    absMin = 12;
    lo = hs + max(round(0.35*stanceLen), absMin);
    hi = hs + round(0.90*stanceLen);
    
    from_idx = max(1, lo - hs);
    if from_idx > numel(heel_vel_z)
        ho = min(max(round(hs + 0.60*stanceLen), lo), hi);
        return
    end
    
    vz_seg = heel_vel_z(from_idx:end);
    thr = mean(vz_seg,'omitnan') + 0.8*std(vz_seg,0,'omitnan');
    
    rel = find(vz_seg > thr, 1, 'first');
    if isempty(rel)
        ho = round(hs + 0.60*stanceLen);
    else
        ho = hs + (from_idx - 1) + rel;
    end
    
    ho = min(max(ho, lo), hi);
end

function mos = calculate_MoS_method1(file_path, hs, ho, to, L0_mm)
    % MÉTHODE 1 : Calcul direct avec filtrage 6 Hz (comme MOS.m)
    
    data = btkReadAcquisition(file_path);
    markers = btkGetMarkers(data);
    marker_names = fieldnames(markers);
    
    % Remplacement M5
    if ~any(strcmp(marker_names, 'RM5')) && any(strcmp(marker_names, 'RM51'))
        markers.RM5 = markers.RM51;
    end
    if ~any(strcmp(marker_names, 'LM5')) && any(strcmp(marker_names, 'LM51'))
        markers.LM5 = markers.LM51;
    end
    
    % Filtrage 6 Hz
    fs = 100; fc = 6;
    [b, a] = butter(2, fc/(fs/2), 'low');
    fns = fieldnames(markers);
    for ii = 1:numel(fns)
        M = markers.(fns{ii});
        if isnumeric(M) && size(M,1) >= 3*max(length(a),length(b))
            markers.(fns{ii}) = filtfilt(b, a, M);
        end
    end
    
    % Détermination du côté
    left_heel_z = markers.LHEE(hs, 3);
    right_heel_z = markers.RHEE(hs, 3);
    
    if left_heel_z < right_heel_z
        side = 'L'; oppside = 'R';
    else
        side = 'R'; oppside = 'L';
    end
    
    % Calcul COM
    COM_x = mean([markers.LPSI(:,1), markers.RPSI(:,1), markers.LASI(:,1), markers.RASI(:,1)], 2);
    COM_y = mean([markers.LPSI(:,2), markers.RPSI(:,2), markers.LASI(:,2), markers.RASI(:,2)], 2);
    COM_z = mean([markers.LPSI(:,3), markers.RPSI(:,3), markers.LASI(:,3), markers.RASI(:,3)], 2);
    
    % Vitesses
    velocities_x = diff(COM_x) * fs;
    velocities_y = diff(COM_y) * fs;
    
    % xCOM
    g = 9810;
    ankle_marker = [side, 'ANK'];
    l_z = abs(markers.(ankle_marker)(1:end-1, 3) - COM_z(1:end-1));
    k = sqrt(g ./ (l_z + 1e-6));
    
    xCOM_x = COM_x(1:end-1) + velocities_x ./ k;
    xCOM_y = COM_y(1:end-1) + velocities_y ./ k;
    
    % Direction ML
    m5_marker_side = [side, 'M5'];
    m5_marker_opp = [oppside, 'M5'];
    
    if (markers.(m5_marker_side)(1, 1) - markers.(m5_marker_opp)(1, 1)) < 0
        directionML = -1;
    else
        directionML = 1;
    end
    
    % Marqueurs
    heel_marker = [side, 'HEE'];
    toe_marker = [side, 'TOE'];
    
    HEEL_y = markers.(heel_marker)(:, 2);
    TOE_y = markers.(toe_marker)(:, 2);
    ANKLE_x = markers.(ankle_marker)(:, 1);
    M5_x = markers.(m5_marker_side)(:, 1);
    
    % Direction AP
    if HEEL_y(to) - HEEL_y(hs) > 0
        directionAP = 1;
    else
        directionAP = -1;
    end
    
    % MoS AP
    MoS_AP_heel = (HEEL_y(hs:ho) - xCOM_y(hs:ho)) * directionAP;
    MoS_AP_toe = (TOE_y(ho:to) - xCOM_y(ho:to)) * directionAP;
    
    % MoS ML
    MoS_ML_ankle = (ANKLE_x(hs:ho) - xCOM_x(hs:ho)) * directionML;
    MoS_ML_m5 = (M5_x(ho:to) - xCOM_x(ho:to)) * directionML;
    mos_ml_total = [MoS_ML_ankle; MoS_ML_m5];
    
    % Résultats
    mos = struct();
    mos.MoS_Heel_Strike_AP = MoS_AP_heel(1);
    mos.MoS_Heel_Strike_ML = mos_ml_total(1);
    mos.MoS_AP_Mean = mean([mean(MoS_AP_heel), mean(MoS_AP_toe)]);
    mos.MoS_ML_Mean = mean(mos_ml_total);
    
    % Normalisation %L0
    if exist('L0_mm','var') && isnumeric(L0_mm) && L0_mm > 0
        mos.MoS_Heel_Strike_AP_P = 100 * (mos.MoS_Heel_Strike_AP / L0_mm);
        mos.MoS_Heel_Strike_ML_P = 100 * (mos.MoS_Heel_Strike_ML / L0_mm);
    end
    
    btkCloseAcquisition(data);
end

function curve = extract_full_mos_curve_method2(base_dir, sujet, surf, essai, hs, ho, to, n_points, L0_mm)
    % MÉTHODE 2 : Avec interpolation (comme Visual_GlobalParticipant_MOS.m)
    
    filename = sprintf('%s_%s_%02d.c3d', sujet, surf, essai);
    c3d_path = fullfile(base_dir, filename);
    
    data = btkReadAcquisition(c3d_path);
    markers = btkGetMarkers(data);
    marker_names = fieldnames(markers);
    
    if ~any(strcmp(marker_names, 'RM5')) && any(strcmp(marker_names, 'RM51'))
        markers.RM5 = markers.RM51;
    end
    if ~any(strcmp(marker_names, 'LM5')) && any(strcmp(marker_names, 'LM51'))
        markers.LM5 = markers.LM51;
    end
    
    % Détermination du côté
    if markers.LHEE(hs,3) < markers.RHEE(hs,3)
        side = 'L'; oppside = 'R';
    else
        side = 'R'; oppside = 'L';
    end
    
    % COM
    COM_AP = mean([markers.LPSI(:,2), markers.RPSI(:,2), markers.LASI(:,2), markers.RASI(:,2)], 2);
    COM_ML = mean([markers.LPSI(:,1), markers.RPSI(:,1), markers.LASI(:,1), markers.RASI(:,1)], 2);
    COM_Z = mean([markers.LPSI(:,3), markers.RPSI(:,3), markers.LASI(:,3), markers.RASI(:,3)], 2);
    
    freqVicon = 100;
    vAP = diff(COM_AP) * freqVicon;
    vML = diff(COM_ML) * freqVicon;
    
    g = 9810;
    ankle_marker = [side, 'ANK'];
    lz = abs(markers.(ankle_marker)(1:end-1, 3) - COM_Z(1:end-1));
    k = sqrt(g ./ (lz + 1e-6));
    
    xCOM_AP = COM_AP(1:end-1) + vAP ./ k;
    xCOM_ML = COM_ML(1:end-1) + vML ./ k;
    
    % Directions
    heel_marker = [side, 'HEE'];
    toe_marker = [side, 'TOE'];
    deltaY = markers.(heel_marker)(to, 2) - markers.(heel_marker)(hs, 2);
    dirAP = sign(deltaY + eps);
    
    m5_side = [side, 'M5']; m5_opp = [oppside, 'M5'];
    dirML = sign(markers.(m5_side)(1,1) - markers.(m5_opp)(1,1) + eps);
    if dirML == 0, dirML = -1; end
    
    % MoS
    HEEL_AP = markers.(heel_marker)(:, 2);
    TOE_AP = markers.(toe_marker)(:, 2);
    ANKLE_ML = markers.(ankle_marker)(:, 1);
    M5_ML = markers.(m5_side)(:, 1);
    
    mos_ap_heel = dirAP * (HEEL_AP(hs:ho) - xCOM_AP(hs:ho));
    mos_ap_toe = dirAP * (TOE_AP(ho:to) - xCOM_AP(ho:to));
    mos_ap_full = [mos_ap_heel; mos_ap_toe];
    
    mos_ml_ankle = (ANKLE_ML(hs:ho) - xCOM_ML(hs:ho)) * dirML;
    mos_ml_m5 = (M5_ML(ho:to) - xCOM_ML(ho:to)) * dirML;
    mos_ml_full = [mos_ml_ankle; mos_ml_m5];
    
    % Interpolation
    x_original = linspace(0, 1, length(mos_ap_full));
    x_norm = linspace(0, 1, n_points);
    
    curve.mos_ap = interp1(x_original, mos_ap_full, x_norm);
    curve.mos_ml = interp1(x_original, mos_ml_full, x_norm);
    
    if exist('L0_mm','var') && isnumeric(L0_mm) && L0_mm > 0
        curve.mos_ap_P = 100 * interp1(x_original, mos_ap_full./L0_mm, x_norm);
        curve.mos_ml_P = 100 * interp1(x_original, mos_ml_full./L0_mm, x_norm);
    else
        curve.mos_ap_P = nan(1, n_points);
        curve.mos_ml_P = nan(1, n_points);
    end
    
    btkCloseAcquisition(data);
end

function mos = calculate_MoS_with_filter(file_path, hs, ho, to, L0_mm, fc, order)
    % Calcul avec paramètres de filtrage personnalisés
    
    data = btkReadAcquisition(file_path);
    markers = btkGetMarkers(data);
    marker_names = fieldnames(markers);
    
    if ~any(strcmp(marker_names, 'RM5')) && any(strcmp(marker_names, 'RM51'))
        markers.RM5 = markers.RM51;
    end
    if ~any(strcmp(marker_names, 'LM5')) && any(strcmp(marker_names, 'LM51'))
        markers.LM5 = markers.LM51;
    end
    
    % Filtrage personnalisé
    fs = 100;
    [b, a] = butter(order, fc/(fs/2), 'low');
    fns = fieldnames(markers);
    for ii = 1:numel(fns)
        M = markers.(fns{ii});
        if isnumeric(M) && size(M,1) >= 3*max(length(a),length(b))
            markers.(fns{ii}) = filtfilt(b, a, M);
        end
    end
    
    % Détermination du côté
    left_heel_z = markers.LHEE(hs, 3);
    right_heel_z = markers.RHEE(hs, 3);
    
    if left_heel_z < right_heel_z
        side = 'L'; oppside = 'R';
    else
        side = 'R'; oppside = 'L';
    end
    
    % Calcul COM
    COM_x = mean([markers.LPSI(:,1), markers.RPSI(:,1), markers.LASI(:,1), markers.RASI(:,1)], 2);
    COM_y = mean([markers.LPSI(:,2), markers.RPSI(:,2), markers.LASI(:,2), markers.RASI(:,2)], 2);
    COM_z = mean([markers.LPSI(:,3), markers.RPSI(:,3), markers.LASI(:,3), markers.RASI(:,3)], 2);
    
    % Vitesses
    velocities_x = diff(COM_x) * fs;
    velocities_y = diff(COM_y) * fs;
    
    % xCOM
    g = 9810;
    ankle_marker = [side, 'ANK'];
    l_z = abs(markers.(ankle_marker)(1:end-1, 3) - COM_z(1:end-1));
    k = sqrt(g ./ (l_z + 1e-6));
    
    xCOM_x = COM_x(1:end-1) + velocities_x ./ k;
    xCOM_y = COM_y(1:end-1) + velocities_y ./ k;
    
    % Direction ML
    m5_marker_side = [side, 'M5'];
    m5_marker_opp = [oppside, 'M5'];
    
    if (markers.(m5_marker_side)(1, 1) - markers.(m5_marker_opp)(1, 1)) < 0
        directionML = -1;
    else
        directionML = 1;
    end
    
    % Marqueurs
    heel_marker = [side, 'HEE'];
    toe_marker = [side, 'TOE'];
    
    HEEL_y = markers.(heel_marker)(:, 2);
    TOE_y = markers.(toe_marker)(:, 2);
    ANKLE_x = markers.(ankle_marker)(:, 1);
    M5_x = markers.(m5_marker_side)(:, 1);
    
    % Direction AP
    if HEEL_y(to) - HEEL_y(hs) > 0
        directionAP = 1;
    else
        directionAP = -1;
    end
    
    % MoS AP
    MoS_AP_heel = (HEEL_y(hs:ho) - xCOM_y(hs:ho)) * directionAP;
    MoS_AP_toe = (TOE_y(ho:to) - xCOM_y(ho:to)) * directionAP;
    
    % MoS ML
    MoS_ML_ankle = (ANKLE_x(hs:ho) - xCOM_x(hs:ho)) * directionML;
    MoS_ML_m5 = (M5_x(ho:to) - xCOM_x(ho:to)) * directionML;
    mos_ml_total = [MoS_ML_ankle; MoS_ML_m5];
    
    % Résultats
    mos = struct();
    mos.MoS_Heel_Strike_AP = MoS_AP_heel(1);
    mos.MoS_Heel_Strike_ML = mos_ml_total(1);
    
    % Normalisation %L0
    if exist('L0_mm','var') && isnumeric(L0_mm) && L0_mm > 0
        mos.MoS_Heel_Strike_AP_P = 100 * (mos.MoS_Heel_Strike_AP / L0_mm);
        mos.MoS_Heel_Strike_ML_P = 100 * (mos.MoS_Heel_Strike_ML / L0_mm);
    end
    
    btkCloseAcquisition(data);
end