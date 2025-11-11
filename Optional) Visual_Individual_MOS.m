%% === VISUALISATION DES MARGES DE STABILITÉ (MoS) ===
% Ce script charge les données MoS, normalise les cycles, et crée
% des graphiques de comparaison pour les 3 surfaces (mm vs %L0)

clear; clc; close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL');
addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\Scripts\btk'));
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));

%% === PARAMÈTRES ===
sujet_id = 'CTL_15';
base_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\jeunes_enfants';
mat_file = sprintf('C:\\Users\\silve\\Desktop\\DOCTORAT\\UNIV MONTREAL\\TRAVAUX-THESE\\Surfaces_Irregulieres\\Datas\\Script\\gaitAnalysisGUI\\result\\MoS\\MoS_results_%s.mat', sujet_id);

% Les 3 surfaces étudiées
surfaces = {'Plat', 'Medium', 'High'};
colors = [0 0.4470 0.7410;    % Bleu pour Plat
          0.4660 0.6740 0.1880; % Vert pour Medium
          0.8500 0.3250 0.0980]; % Rouge pour High

% Nombre de points pour la normalisation
n_points = 100;

%% === L0 du participant (pour %L0) ===
% Map attendue: containers.Map key=participant (ex 'CTL_01'), value=L0 en mètres
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL')
load('l0_participants.mat', 'l0_map');
if ~isKey(l0_map, sujet_id)
    error('Participant %s non trouvé dans l0_map.', sujet_id);
end
L0_mm = l0_map(sujet_id) * 1000;  % conversion m -> mm
fprintf('Participant %s: L0 = %.1f mm\n\n', sujet_id, L0_mm);

%% === CHARGEMENT DES DONNÉES ===
fprintf('📂 Chargement des données MoS pour %s...\n', sujet_id);

if ~isfile(mat_file)
    error('❌ Fichier MAT introuvable : %s', mat_file);
end

load(mat_file, 'MoS_data');
results = MoS_data.results;

fprintf('✅ Données chargées : %d cycles au total\n\n', height(results));

%% === EXTRACTION ET NORMALISATION DES COURBES MoS ===
fprintf('🔄 Extraction des courbes MoS complètes...\n');

% Initialisation des structures pour stocker les courbes par surface
mos_ap   = cell(3, 1);     % en mm
mos_ml   = cell(3, 1);     % en mm
mos_ap_P = cell(3, 1);     % en %L0
mos_ml_P = cell(3, 1);     % en %L0

for s = 1:length(surfaces)
    surf = surfaces{s};
    
    % Filtrer par surface
    idx_surf = strcmp(results.Surface, surf);
    results_surf = results(idx_surf, :);
    
    fprintf('Surface %s : %d cycles\n', surf, height(results_surf));
    
    % Extraction des courbes
    mos_ap_temp   = [];
    mos_ml_temp   = [];
    mos_apP_temp  = [];
    mos_mlP_temp  = [];
    
    for i = 1:height(results_surf)
        row = results_surf(i, :);
        try
            curve = extract_full_mos_curve(base_dir, row, n_points, L0_mm);
            mos_ap_temp  = [mos_ap_temp;  curve.mos_ap];
            mos_ml_temp  = [mos_ml_temp;  curve.mos_ml];
            mos_apP_temp = [mos_apP_temp; curve.mos_ap_P];
            mos_mlP_temp = [mos_mlP_temp; curve.mos_ml_P];
        catch ME
            fprintf('   ⚠️  Erreur cycle %d: %s\n', i, ME.message);
        end
    end
    
    mos_ap{s}    = mos_ap_temp;
    mos_ml{s}    = mos_ml_temp;
    mos_ap_P{s}  = mos_apP_temp;
    mos_ml_P{s}  = mos_mlP_temp;
    
    fprintf('   ✓ %d courbes extraites\n', size(mos_ap_temp, 1));
end

fprintf('\n✅ Extraction terminée\n\n');

% axe normalisé (0-100 %) → on en aura besoin pour la structure
x_norm = linspace(0, 100, n_points);

%% === STRUCTURE COMPLÈTE DES CYCLES MoS (AP & ML, mm & %L0, toutes surfaces) ===
MoS_store = struct();
MoS_store.Sujet     = sujet_id;
MoS_store.L0_mm     = L0_mm;
MoS_store.n_points  = n_points;
MoS_store.x_norm    = x_norm;
MoS_store.Surfaces  = surfaces;

for s = 1:length(surfaces)
    surf = surfaces{s};

    % on crée le sous-struct surface
    MoS_store.(surf) = struct();

    % ===================== AP =====================
    % --- AP en mm ---
    all_ap_mm = mos_ap{s};
    MoS_store.(surf).AP.mm.all       = all_ap_mm;
    MoS_store.(surf).AP.mm.n_cycles  = size(all_ap_mm, 1);
    if ~isempty(all_ap_mm)
        MoS_store.(surf).AP.mm.mean  = mean(all_ap_mm, 1);
        MoS_store.(surf).AP.mm.std   = std (all_ap_mm, 0, 1);
    else
        MoS_store.(surf).AP.mm.mean  = [];
        MoS_store.(surf).AP.mm.std   = [];
    end

    % --- AP en %L0 ---
    all_ap_pL0 = mos_ap_P{s};
    MoS_store.(surf).AP.pL0.all      = all_ap_pL0;
    MoS_store.(surf).AP.pL0.n_cycles = size(all_ap_pL0, 1);
    if ~isempty(all_ap_pL0)
        MoS_store.(surf).AP.pL0.mean = mean(all_ap_pL0, 1);
        MoS_store.(surf).AP.pL0.std  = std (all_ap_pL0, 0, 1);
    else
        MoS_store.(surf).AP.pL0.mean = [];
        MoS_store.(surf).AP.pL0.std  = [];
    end

    % ===================== ML =====================
    % --- ML en mm ---
    all_ml_mm = mos_ml{s};
    MoS_store.(surf).ML.mm.all       = all_ml_mm;
    MoS_store.(surf).ML.mm.n_cycles  = size(all_ml_mm, 1);
    if ~isempty(all_ml_mm)
        MoS_store.(surf).ML.mm.mean  = mean(all_ml_mm, 1);
        MoS_store.(surf).ML.mm.std   = std (all_ml_mm, 0, 1);
    else
        MoS_store.(surf).ML.mm.mean  = [];
        MoS_store.(surf).ML.mm.std   = [];
    end

    % --- ML en %L0 ---
    all_ml_pL0 = mos_ml_P{s};
    MoS_store.(surf).ML.pL0.all      = all_ml_pL0;
    MoS_store.(surf).ML.pL0.n_cycles = size(all_ml_pL0, 1);
    if ~isempty(all_ml_pL0)
        MoS_store.(surf).ML.pL0.mean = mean(all_ml_pL0, 1);
        MoS_store.(surf).ML.pL0.std  = std (all_ml_pL0, 0, 1);
    else
        MoS_store.(surf).ML.pL0.mean = [];
        MoS_store.(surf).ML.pL0.std  = [];
    end
end

% === SAUVEGARDE ===
save_dir_struct = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\MoS\Cycles';
if ~exist(save_dir_struct, 'dir')
    mkdir(save_dir_struct);
end

save(fullfile(save_dir_struct, sprintf('MoS_store_%s.mat', sujet_id)), 'MoS_store');
fprintf('✅ Structure complète MoS sauvegardée pour %s\n', sujet_id);

%% === VISUALISATION MoS ANTÉRO-POSTÉRIEUR (mm vs %L0) ===
fprintf('📊 Création du graphique MoS AP (mm et %%L0)...\n');

figure('Position', [100, 100, 1400, 550]);
tiledlayout(1,2, 'Padding','compact','TileSpacing','compact');

% --- (A) AP en mm
nexttile; hold on;
for s = 1:length(surfaces)
    if ~isempty(mos_ap{s})
        mean_ap = mean(mos_ap{s}, 1);
        std_ap  = std (mos_ap{s}, 0, 1);
        plot(x_norm, mean_ap, '-', 'LineWidth', 2.2, ...
             'Color', colors(s,:), 'DisplayName', surfaces{s});
        fill([x_norm, fliplr(x_norm)], ...
             [mean_ap + std_ap, fliplr(mean_ap - std_ap)], ...
             colors(s,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
    end
end
xlabel('% du cycle d''appui'); ylabel('MoS AP (mm)');
title(sprintf('MoS AP (mm) – %s', sujet_id));
legend('Location','best'); grid on; box on; hold off;

% --- (B) AP en %L0
nexttile; hold on;
for s = 1:length(surfaces)
    if ~isempty(mos_ap_P{s})
        mean_apP = mean(mos_ap_P{s}, 1);
        std_apP  = std (mos_ap_P{s}, 0, 1);
        plot(x_norm, mean_apP, '--', 'LineWidth', 2.2, ...
             'Color', colors(s,:), 'DisplayName', [surfaces{s} ' (%L0)']);
        fill([x_norm, fliplr(x_norm)], ...
             [mean_apP + std_apP, fliplr(mean_apP - std_apP)], ...
             colors(s,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
    end
end
xlabel('% du cycle d''appui'); ylabel('MoS AP (%L0)');
title(sprintf('MoS AP normalisée (%%L0) – %s', sujet_id));
legend('Location','best'); grid on; box on; hold off;

save_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Fig\MOS';
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
saveas(gcf, fullfile(save_dir, sprintf('MoS_AP_%s.png', sujet_id)));
fprintf('✅ Figure MoS AP sauvegardée : %s\n', save_dir);

%% === VISUALISATION MoS MÉDIO-LATÉRAL (mm vs %L0) ===
fprintf('📊 Création du graphique MoS ML (mm et %%L0)...\n');

figure('Position', [150, 150, 1400, 550]);
tiledlayout(1,2, 'Padding','compact','TileSpacing','compact');

% --- (A) ML en mm
nexttile; hold on;
for s = 1:length(surfaces)
    if ~isempty(mos_ml{s})
        mean_ml = mean(mos_ml{s}, 1);
        std_ml  = std (mos_ml{s}, 0, 1);
        plot(x_norm, mean_ml, '-', 'LineWidth', 2.2, ...
             'Color', colors(s,:), 'DisplayName', surfaces{s});
        fill([x_norm, fliplr(x_norm)], ...
             [mean_ml + std_ml, fliplr(mean_ml - std_ml)], ...
             colors(s,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
    end
end
xlabel('% du cycle d''appui'); ylabel('MoS ML (mm)');
title(sprintf('MoS ML (mm) – %s', sujet_id));
legend('Location','best'); grid on; box on; hold off;

% --- (B) ML en %L0
nexttile; hold on;
for s = 1:length(surfaces)
    if ~isempty(mos_ml_P{s})
        mean_mlP = mean(mos_ml_P{s}, 1);
        std_mlP  = std (mos_ml_P{s}, 0, 1);
        plot(x_norm, mean_mlP, '--', 'LineWidth', 2.2, ...
             'Color', colors(s,:), 'DisplayName', [surfaces{s} ' (%L0)']);
        fill([x_norm, fliplr(x_norm)], ...
             [mean_mlP + std_mlP, fliplr(mean_mlP - std_mlP)], ...
             colors(s,:), 'FaceAlpha', 0.18, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
    end
end
xlabel('% du cycle d''appui'); ylabel('MoS ML (%L0)');
title(sprintf('MoS ML normalisée (%%L0) – %s', sujet_id));
legend('Location','best'); grid on; box on; hold off;

% === ENREGISTREMENT ===
saveas(gcf, fullfile(save_dir, sprintf('MoS_ML_%s.png', sujet_id)));
fprintf('✅ Figure MoS ML sauvegardée : %s\n', save_dir);

%% === STATISTIQUES DESCRIPTIVES (brut et %L0) ===
fprintf('\n📊 STATISTIQUES DESCRIPTIVES (brut et %%L0)\n');
fprintf('══════════════════════════════════════════════════════════\n');

for s = 1:length(surfaces)
    fprintf('\n%s :\n', surfaces{s});
    if ~isempty(mos_ap{s})
        fprintf('  MoS AP (mm)  : %.2f ± %.2f | min %.2f | max %.2f\n', ...
            mean(mos_ap{s}(:)), std(mos_ap{s}(:)), min(mos_ap{s}(:)), max(mos_ap{s}(:)));
    end
    if ~isempty(mos_ap_P{s})
        fprintf('  MoS AP (%%L0) : %.2f ± %.2f | min %.2f | max %.2f\n', ...
            mean(mos_ap_P{s}(:)), std(mos_ap_P{s}(:)), min(mos_ap_P{s}(:)), max(mos_ap_P{s}(:)));
    end
    if ~isempty(mos_ml{s})
        fprintf('  MoS ML (mm)  : %.2f ± %.2f | min %.2f | max %.2f\n', ...
            mean(mos_ml{s}(:)), std(mos_ml{s}(:)), min(mos_ml{s}(:)), max(mos_ml{s}(:)));
    end
    if ~isempty(mos_ml_P{s})
        fprintf('  MoS ML (%%L0) : %.2f ± %.2f | min %.2f | max %.2f\n', ...
            mean(mos_ml_P{s}(:)), std(mos_ml_P{s}(:)), min(mos_ml_P{s}(:)), max(mos_ml_P{s}(:)));
    end
end

fprintf('\n✅ Visualisations créées avec succès!\n');

%% ========== FONCTION ==========

function curve = extract_full_mos_curve(base_dir, row, n_points, L0_mm)
    % Extrait et normalise la courbe complète de MoS (AP et ML) pour un cycle
    % Renvoie aussi la version normalisée en %L0 (si L0_mm > 0)
    % Repère Vicon:
    %   Y = AP (direction de marche), X = ML, Z = vertical

    % Fichier C3D
    filename = sprintf('%s_%s_%02d.c3d', row.Sujet{1}, row.Surface{1}, row.Essai);
    c3d_path = fullfile(base_dir, filename);

    % Événements
    hs = row.HS_Frame; ho = row.HO_Frame; to = row.TO_Frame;

    % Lecture C3D
    data = btkReadAcquisition(c3d_path);
    markers = btkGetMarkers(data);
    marker_names = fieldnames(markers);

    % Gestion M5/M51
    if ~any(strcmp(marker_names, 'RM5')) && any(strcmp(marker_names, 'RM51')), markers.RM5 = markers.RM51; end
    if ~any(strcmp(marker_names, 'LM5')) && any(strcmp(marker_names, 'LM51')), markers.LM5 = markers.LM51; end

    % Côté
    if strcmp(row.Cote{1}, 'Gauche'), side = 'L'; oppside = 'R'; else, side = 'R'; oppside = 'L'; end

    % COM (AP=Y, ML=X)
    COM_AP = mean([markers.LPSI(:,2), markers.RPSI(:,2), markers.LASI(:,2), markers.RASI(:,2)], 2);
    COM_ML = mean([markers.LPSI(:,1), markers.RPSI(:,1), markers.LASI(:,1), markers.RASI(:,1)], 2);
    COM_Z  = mean([markers.LPSI(:,3), markers.RPSI(:,3), markers.LASI(:,3), markers.RASI(:,3)], 2);

    % Vitesses
    freqVicon = 100;
    vAP = diff(COM_AP) * freqVicon;
    vML = diff(COM_ML) * freqVicon;

    % xCOM
    g = 9810;
    ankle_marker = [side, 'ANK'];
    lz = abs(markers.(ankle_marker)(1:end-1, 3) - COM_Z(1:end-1));
    k  = sqrt(g ./ (lz + 1e-6));
    xCOM_AP = COM_AP(1:end-1) + vAP ./ k;
    xCOM_ML = COM_ML(1:end-1) + vML ./ k;

    % Sens AP (optionnel/stable)
    heel_marker = [side, 'HEE'];
    deltaY = markers.(heel_marker)(to, 2) - markers.(heel_marker)(hs, 2);
    if deltaY > 0, dirAP = 1; else, dirAP = -1; end

    % Direction ML
    m5_side = [side, 'M5']; m5_opp = [oppside, 'M5'];
    if (markers.(m5_side)(1,1) - markers.(m5_opp)(1,1)) < 0, dirML = -1; else, dirML = 1; end

    % Marqueurs utiles
    toe_marker = [side, 'TOE'];
    HEEL_AP  = markers.(heel_marker)(:, 2);
    TOE_AP   = markers.(toe_marker)(:, 2);
    ANKLE_ML = markers.(ankle_marker)(:, 1);
    M5_ML    = markers.(m5_side)(:, 1);

    % Courbes MoS (mm)
    mos_ap_heel = dirAP * (HEEL_AP(hs:ho) - xCOM_AP(hs:ho));
    mos_ap_toe  = dirAP * (TOE_AP(ho:to)  - xCOM_AP(ho:to));
    mos_ap_full = [mos_ap_heel; mos_ap_toe];

    mos_ml_ankle = (ANKLE_ML(hs:ho) - xCOM_ML(hs:ho)) * dirML;
    mos_ml_m5    = (M5_ML(ho:to)   - xCOM_ML(ho:to))  * dirML;
    mos_ml_full  = [mos_ml_ankle; mos_ml_m5];

    % Interpolation sur 0–100 %
    x_original = linspace(0, 1, length(mos_ap_full));
    x_norm = linspace(0, 1, n_points);
    curve.mos_ap = interp1(x_original, mos_ap_full, x_norm);
    curve.mos_ml = interp1(x_original, mos_ml_full, x_norm);

    % Séries normalisées en %L0
    if exist('L0_mm','var') && isnumeric(L0_mm) && isfinite(L0_mm) && L0_mm > 0
        curve.mos_ap_P = 100 * interp1(x_original, mos_ap_full./L0_mm, x_norm);
        curve.mos_ml_P = 100 * interp1(x_original, mos_ml_full ./L0_mm, x_norm);
    else
        % Si L0 absent/invalide, renvoie des NaN pour garder la compatibilité
        curve.mos_ap_P = nan(1, n_points);
        curve.mos_ml_P = nan(1, n_points);
        warning('L0_mm invalide ou manquant: les séries %%L0 sont NaN.');
    end

    btkCloseAcquisition(data);
end