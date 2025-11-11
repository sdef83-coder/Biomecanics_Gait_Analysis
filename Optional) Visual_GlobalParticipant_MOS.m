%% ================================================================
%  VISUALISATION / AGRÉGATION DES MARGES DE STABILITÉ (MoS)
%  À PARTIR DES FICHIERS MoS_results_XXX.mat + C3D
%  - Même logique que ton script individuel
%  - Mais appliqué à tous les participants de chaque groupe
%  - Sort 8 figures :
%     (1) Effet âge – AP mm
%     (2) Effet âge – AP %L0
%     (3) Effet âge – ML mm
%     (4) Effet âge – ML %L0
%     (5) Effet surface – AP mm
%     (6) Effet surface – AP %L0
%     (7) Effet surface – ML mm
%     (8) Effet surface – ML %L0
% ================================================================
clear; clc; close all;

%% === PATHS GÉNÉRAUX ===
root_res_all = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL';
cd(root_res_all);

addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\Scripts\btk'));
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI'));

% dossier où sont les MoS_results_XX.mat
mos_res_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\MoS';

% dossier de sauvegarde des figures
save_fig_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Fig\MOS_GROUPES';
if ~exist(save_fig_dir, 'dir'); mkdir(save_fig_dir); end

%% === GROUPES DE PARTICIPANTS ===
ParticipantGroup
group_names = fieldnames(Group);

%% === PARAMÈTRES COMMUNS ===
surfaces = {'Plat', 'Medium', 'High'};
n_points = 100;             % normalisation 0-100 %
x_norm   = linspace(0,100,n_points);

% couleurs demandées
color_groups = [
    0    0.6   1;    % JeunesEnfants
    1    0.5   0;    % Enfants
    0    0.7   0.3;  % Adolescents
    0.5  0     1];   % Adultes

color_surfaces = [
    0.2 0.4 0.8;  % Bleu - Plat
    0.2 0.7 0.3;  % Vert - Medium
    0.8 0.2 0.2]; % High   -> rouge

%% === L0 DES PARTICIPANTS ===
% même logique que ton script de base
load('l0_participants.mat', 'l0_map');

%% === STRUCTURE POUR STOCKER LES MOYENNES DE GROUPE ===
% groupData.(Groupe).(Surface).AP_mm_mean, etc.
groupData = struct();
groupData.x_norm = x_norm;

fprintf('🔄 DÉBUT DU TRAITEMENT PAR GROUPE...\n');

for g = 1:numel(group_names)
    gname = group_names{g};
    sujets = Group.(gname);
    fprintf('\n===== GROUPE : %s =====\n', gname);

    % choisir le dossier des C3D selon le groupe
    switch gname
        case 'JeunesEnfants'
            base_dir_group = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\jeunes_enfants';
        case 'Enfants'
            base_dir_group = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\enfants';
        case 'Adolescents'
            base_dir_group = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\adolescents';
        case 'Adultes'
            base_dir_group = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\adults';
        otherwise
            error('Dossier de données inconnu pour le groupe %s', gname);
    end

    % initialiser les contenants pour empiler "1 courbe moyenne PAR SUJET"
    for s = 1:numel(surfaces)
        surf = surfaces{s};
        groupData.(gname).(surf).AP_mm  = [];  % chaque ligne = 1 sujet
        groupData.(gname).(surf).AP_pL0 = [];
        groupData.(gname).(surf).ML_mm  = [];
        groupData.(gname).(surf).ML_pL0 = [];
    end

    % === boucle sujets du groupe ===
    for iS = 1:numel(sujets)
        sujet_id = sujets{iS};
        fprintf(' → Sujet %s\n', sujet_id);

        % récupérer L0
        if ~isKey(l0_map, sujet_id)
            fprintf('   ⚠️ L0 manquant pour %s, sujet ignoré\n', sujet_id);
            continue;
        end
        L0_mm = l0_map(sujet_id) * 1000;

        % charger le .mat de MoS
        mat_file = fullfile(mos_res_dir, sprintf('MoS_results_%s.mat', sujet_id));
        if ~isfile(mat_file)
            fprintf('   ⚠️ Fichier MoS introuvable pour %s, ignoré\n', sujet_id);
            continue;
        end
        tmp = load(mat_file, 'MoS_data');
        results = tmp.MoS_data.results;

        % pour stocker les courbes de ce SUJET (moyennées ensuite)
        subj_AP_mm  = struct(); subj_AP_pL0 = struct();
        subj_ML_mm  = struct(); subj_ML_pL0 = struct();

        for s = 1:numel(surfaces)
            surf = surfaces{s};
            idx_surf = strcmp(results.Surface, surf);
            res_surf = results(idx_surf, :);

            mos_ap_temp  = [];
            mos_ml_temp  = [];
            mos_apP_temp = [];
            mos_mlP_temp = [];

            for i = 1:height(res_surf)
                row = res_surf(i, :);
                try
                    curve = extract_full_mos_curve(base_dir_group, row, n_points, L0_mm);
                    mos_ap_temp  = [mos_ap_temp;  curve.mos_ap];
                    mos_ml_temp  = [mos_ml_temp;  curve.mos_ml];
                    mos_apP_temp = [mos_apP_temp; curve.mos_ap_P];
                    mos_mlP_temp = [mos_mlP_temp; curve.mos_ml_P];
                catch ME
                    fprintf('     ⚠️ Erreur cycle %d (%s, %s): %s\n', i, sujet_id, surf, ME.message);
                end
            end

            % si pas de cycle pour cette surface → on ne met rien
            if isempty(mos_ap_temp)
                continue;
            end

            % → moyenne PAR SUJET (ce que tu voulais pour ne pas surreprésenter)
            subj_AP_mm.(surf)  = mean(mos_ap_temp, 1);
            subj_AP_pL0.(surf) = mean(mos_apP_temp, 1);
            subj_ML_mm.(surf)  = mean(mos_ml_temp, 1);
            subj_ML_pL0.(surf) = mean(mos_mlP_temp, 1);

            % → on empile dans le groupe
            groupData.(gname).(surf).AP_mm  = [groupData.(gname).(surf).AP_mm;  subj_AP_mm.(surf)];
            groupData.(gname).(surf).AP_pL0 = [groupData.(gname).(surf).AP_pL0; subj_AP_pL0.(surf)];
            groupData.(gname).(surf).ML_mm  = [groupData.(gname).(surf).ML_mm;  subj_ML_mm.(surf)];
            groupData.(gname).(surf).ML_pL0 = [groupData.(gname).(surf).ML_pL0; subj_ML_pL0.(surf)];
        end
    end

    % === moyenne DE GROUPE une fois tous les sujets passés ===
    for s = 1:numel(surfaces)
        surf = surfaces{s};

        A = groupData.(gname).(surf).AP_mm;
        if ~isempty(A)
            groupData.(gname).(surf).AP_mm_mean = mean(A, 1);
            groupData.(gname).(surf).AP_mm_std  = std (A, 0, 1);
        else
            groupData.(gname).(surf).AP_mm_mean = [];
            groupData.(gname).(surf).AP_mm_std  = [];
        end

        A = groupData.(gname).(surf).AP_pL0;
        if ~isempty(A)
            groupData.(gname).(surf).AP_pL0_mean = mean(A, 1);
            groupData.(gname).(surf).AP_pL0_std  = std (A, 0, 1);
        else
            groupData.(gname).(surf).AP_pL0_mean = [];
            groupData.(gname).(surf).AP_pL0_std  = [];
        end

        A = groupData.(gname).(surf).ML_mm;
        if ~isempty(A)
            groupData.(gname).(surf).ML_mm_mean = mean(A, 1);
            groupData.(gname).(surf).ML_mm_std  = std (A, 0, 1);
        else
            groupData.(gname).(surf).ML_mm_mean = [];
            groupData.(gname).(surf).ML_mm_std  = [];
        end

        A = groupData.(gname).(surf).ML_pL0;
        if ~isempty(A)
            groupData.(gname).(surf).ML_pL0_mean = mean(A, 1);
            groupData.(gname).(surf).ML_pL0_std  = std (A, 0, 1);
        else
            groupData.(gname).(surf).ML_pL0_mean = [];
            groupData.(gname).(surf).ML_pL0_std  = [];
        end
    end
end

fprintf('\n✅ AGRÉGATION FINIE.\n');

%% === SAUVEGARDE DES COURBES DE GROUPE ===
% On sauvegarde tout ce qu'on vient de calculer dans un .mat
% pour pouvoir le recharger dans un autre script (ANOVA, stats, export, etc.)

save_group_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\MoS\Groupes';
if ~exist(save_group_dir, 'dir'); mkdir(save_group_dir); end

groupMoS = struct();
groupMoS.groupData  = groupData;      % tout ce qu'on a agrégé
groupMoS.surfaces   = surfaces;
groupMoS.groupNames = group_names;
groupMoS.n_points   = n_points;
groupMoS.x_norm     = x_norm;
groupMoS.colors.groups   = color_groups;
groupMoS.colors.surfaces = color_surfaces;

save(fullfile(save_group_dir, 'MoS_groupes.mat'), 'groupMoS');

fprintf('💾 Données de groupes sauvegardées dans : %s\\MoS_groupes.mat\n', save_group_dir);

%% ================================================================
%  FIGURES 1–4 : EFFET ÂGE
%  → 3 subplots : Plat / Medium / High
%  avec bandes ±1 SD
% ================================================================
showSTD = true;   % mets false si tu veux enlever les bandes

% ----- (1) Effet âge – AP (mm) -----
f1 = figure('Name','Effet age - MoS AP (mm)','Position',[100 100 1300 550]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for s = 1:numel(surfaces)
    nexttile; hold on;
    for g = 1:numel(group_names)
        gname  = group_names{g};
        y_mean = groupData.(gname).(surfaces{s}).AP_mm_mean;
        y_std  = groupData.(gname).(surfaces{s}).AP_mm_std;

        if isempty(y_mean), continue; end

        % bande ± SD
        if showSTD && ~isempty(y_std)
            fill([x_norm, fliplr(x_norm)], ...
                 [y_mean + y_std, fliplr(y_mean - y_std)], ...
                 color_groups(g,:), ...
                 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end

        % courbe moyenne
        plot(x_norm, y_mean, 'LineWidth', 2.1, ...
            'Color', color_groups(g,:), ...
            'DisplayName', gname);
    end
    title(sprintf('Surface : %s', surfaces{s}));
    xlabel('% du cycle d''appui'); ylabel('MoS AP (mm)');
    grid on; box on;
    if s == 1, legend('Location','best'); end
end
saveas(f1, fullfile(save_fig_dir, 'EffetAge_MoS_AP_mm.png'));

% ----- (2) Effet âge – AP (%L0) -----
f2 = figure('Name','Effet age - MoS AP (%L0)','Position',[120 120 1300 550]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for s = 1:numel(surfaces)
    nexttile; hold on;
    for g = 1:numel(group_names)
        gname  = group_names{g};
        y_mean = groupData.(gname).(surfaces{s}).AP_pL0_mean;
        y_std  = groupData.(gname).(surfaces{s}).AP_pL0_std;

        if isempty(y_mean), continue; end

        if showSTD && ~isempty(y_std)
            fill([x_norm, fliplr(x_norm)], ...
                 [y_mean + y_std, fliplr(y_mean - y_std)], ...
                 color_groups(g,:), ...
                 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end

        plot(x_norm, y_mean, '--', 'LineWidth', 2.1, ...
            'Color', color_groups(g,:), ...
            'DisplayName', gname);
    end
    title(sprintf('Surface : %s', surfaces{s}));
    xlabel('% du cycle d''appui'); ylabel('MoS AP (%L0)');
    grid on; box on;
    if s == 1, legend('Location','best'); end
end
saveas(f2, fullfile(save_fig_dir, 'EffetAge_MoS_AP_pL0.png'));

% ----- (3) Effet âge – ML (mm) -----
f3 = figure('Name','Effet age - MoS ML (mm)','Position',[140 140 1300 550]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for s = 1:numel(surfaces)
    nexttile; hold on;
    for g = 1:numel(group_names)
        gname  = group_names{g};
        y_mean = groupData.(gname).(surfaces{s}).ML_mm_mean;
        y_std  = groupData.(gname).(surfaces{s}).ML_mm_std;

        if isempty(y_mean), continue; end

        if showSTD && ~isempty(y_std)
            fill([x_norm, fliplr(x_norm)], ...
                 [y_mean + y_std, fliplr(y_mean - y_std)], ...
                 color_groups(g,:), ...
                 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end

        plot(x_norm, y_mean, 'LineWidth', 2.1, ...
            'Color', color_groups(g,:), ...
            'DisplayName', gname);
    end
    title(sprintf('Surface : %s', surfaces{s}));
    xlabel('% du cycle d''appui'); ylabel('MoS ML (mm)');
    grid on; box on;
    if s == 1, legend('Location','best'); end
end
saveas(f3, fullfile(save_fig_dir, 'EffetAge_MoS_ML_mm.png'));

% ----- (4) Effet âge – ML (%L0) -----
f4 = figure('Name','Effet age - MoS ML (%L0)','Position',[160 160 1300 550]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for s = 1:numel(surfaces)
    nexttile; hold on;
    for g = 1:numel(group_names)
        gname  = group_names{g};
        y_mean = groupData.(gname).(surfaces{s}).ML_pL0_mean;
        y_std  = groupData.(gname).(surfaces{s}).ML_pL0_std;

        if isempty(y_mean), continue; end

        if showSTD && ~isempty(y_std)
            fill([x_norm, fliplr(x_norm)], ...
                 [y_mean + y_std, fliplr(y_mean - y_std)], ...
                 color_groups(g,:), ...
                 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end

        plot(x_norm, y_mean, '--', 'LineWidth', 2.1, ...
            'Color', color_groups(g,:), ...
            'DisplayName', gname);
    end
    title(sprintf('Surface : %s', surfaces{s}));
    xlabel('% du cycle d''appui'); ylabel('MoS ML (%L0)');
    grid on; box on;
    if s == 1, legend('Location','best'); end
end
saveas(f4, fullfile(save_fig_dir, 'EffetAge_MoS_ML_pL0.png'));

% ================================================================
%  FIGURES 5–8 : EFFET SURFACE
%  → 4 subplots : 1 par groupe
%  avec bandes ±1 SD
% ================================================================

% ----- (5) Effet surface – AP (mm) -----
f5 = figure('Name','Effet surface - MoS AP (mm)','Position',[180 180 1400 650]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for g = 1:numel(group_names)
    gname = group_names{g};
    nexttile; hold on;

    for s = 1:numel(surfaces)
        y_mean = groupData.(gname).(surfaces{s}).AP_mm_mean;
        y_std  = groupData.(gname).(surfaces{s}).AP_mm_std;

        if isempty(y_mean), continue; end

        if showSTD && ~isempty(y_std)
            fill([x_norm, fliplr(x_norm)], ...
                 [y_mean + y_std, fliplr(y_mean - y_std)], ...
                 color_surfaces(s,:), ...
                 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end

        plot(x_norm, y_mean, 'LineWidth', 2.1, ...
            'Color', color_surfaces(s,:), ...
            'DisplayName', surfaces{s});
    end

    title(gname);
    xlabel('% du cycle d''appui'); ylabel('MoS AP (mm)');
    grid on; box on;
    if g == 1, legend('Location','best'); end
end
saveas(f5, fullfile(save_fig_dir, 'EffetSurface_MoS_AP_mm.png'));

% ----- (6) Effet surface – AP (%L0) -----
f6 = figure('Name','Effet surface - MoS AP (%L0)','Position',[200 200 1400 650]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for g = 1:numel(group_names)
    gname = group_names{g};
    nexttile; hold on;

    for s = 1:numel(surfaces)
        y_mean = groupData.(gname).(surfaces{s}).AP_pL0_mean;
        y_std  = groupData.(gname).(surfaces{s}).AP_pL0_std;

        if isempty(y_mean), continue; end

        if showSTD && ~isempty(y_std)
            fill([x_norm, fliplr(x_norm)], ...
                 [y_mean + y_std, fliplr(y_mean - y_std)], ...
                 color_surfaces(s,:), ...
                 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end

        plot(x_norm, y_mean, '--', 'LineWidth', 2.1, ...
            'Color', color_surfaces(s,:), ...
            'DisplayName', surfaces{s});
    end

    title(gname);
    xlabel('% du cycle d''appui'); ylabel('MoS AP (%L0)');
    grid on; box on;
    if g == 1, legend('Location','best'); end
end
saveas(f6, fullfile(save_fig_dir, 'EffetSurface_MoS_AP_pL0.png'));

% ----- (7) Effet surface – ML (mm) -----
f7 = figure('Name','Effet surface - MoS ML (mm)','Position',[220 220 1400 650]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for g = 1:numel(group_names)
    gname = group_names{g};
    nexttile; hold on;

    for s = 1:numel(surfaces)
        y_mean = groupData.(gname).(surfaces{s}).ML_mm_mean;
        y_std  = groupData.(gname).(surfaces{s}).ML_mm_std;

        if isempty(y_mean), continue; end

        if showSTD && ~isempty(y_std)
            fill([x_norm, fliplr(x_norm)], ...
                 [y_mean + y_std, fliplr(y_mean - y_std)], ...
                 color_surfaces(s,:), ...
                 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end

        plot(x_norm, y_mean, 'LineWidth', 2.1, ...
            'Color', color_surfaces(s,:), ...
            'DisplayName', surfaces{s});
    end

    title(gname);
    xlabel('% du cycle d''appui'); ylabel('MoS ML (mm)');
    grid on; box on;
    if g == 1, legend('Location','best'); end
end
saveas(f7, fullfile(save_fig_dir, 'EffetSurface_MoS_ML_mm.png'));

% ----- (8) Effet surface – ML (%L0) -----
f8 = figure('Name','Effet surface - MoS ML (%L0)','Position',[240 240 1400 650]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

for g = 1:numel(group_names)
    gname = group_names{g};
    nexttile; hold on;

    for s = 1:numel(surfaces)
        y_mean = groupData.(gname).(surfaces{s}).ML_pL0_mean;
        y_std  = groupData.(gname).(surfaces{s}).ML_pL0_std;

        if isempty(y_mean), continue; end

        if showSTD && ~isempty(y_std)
            fill([x_norm, fliplr(x_norm)], ...
                 [y_mean + y_std, fliplr(y_mean - y_std)], ...
                 color_surfaces(s,:), ...
                 'FaceAlpha', 0.13, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end

        plot(x_norm, y_mean, '--', 'LineWidth', 2.1, ...
            'Color', color_surfaces(s,:), ...
            'DisplayName', surfaces{s});
    end

    title(gname);
    xlabel('% du cycle d''appui'); ylabel('MoS ML (%L0)');
    grid on; box on;
    if g == 1, legend('Location','best'); end
end
saveas(f8, fullfile(save_fig_dir, 'EffetSurface_MoS_ML_pL0.png'));

fprintf('\n🎯 8 figures générées dans : %s\n', save_fig_dir);


%% ========== FONCTION (TA MÊME LOGIQUE) ==========
function curve = extract_full_mos_curve(base_dir, row, n_points, L0_mm)
    % même que ton script de base

    filename = sprintf('%s_%s_%02d.c3d', row.Sujet{1}, row.Surface{1}, row.Essai);
    c3d_path = fullfile(base_dir, filename);

    hs = row.HS_Frame; ho = row.HO_Frame; to = row.TO_Frame;

    data = btkReadAcquisition(c3d_path);
    markers = btkGetMarkers(data);
    marker_names = fieldnames(markers);

    if ~any(strcmp(marker_names, 'RM5')) && any(strcmp(marker_names, 'RM51')), markers.RM5 = markers.RM51; end
    if ~any(strcmp(marker_names, 'LM5')) && any(strcmp(marker_names, 'LM51')), markers.LM5 = markers.LM51; end

    if strcmp(row.Cote{1}, 'Gauche'), side = 'L'; oppside = 'R'; else, side = 'R'; oppside = 'L'; end

    COM_AP = mean([markers.LPSI(:,2), markers.RPSI(:,2), markers.LASI(:,2), markers.RASI(:,2)], 2);
    COM_ML = mean([markers.LPSI(:,1), markers.RPSI(:,1), markers.LASI(:,1), markers.RASI(:,1)], 2);
    COM_Z  = mean([markers.LPSI(:,3), markers.RPSI(:,3), markers.LASI(:,3), markers.RASI(:,3)], 2);

    freqVicon = 100;
    vAP = diff(COM_AP) * freqVicon;
    vML = diff(COM_ML) * freqVicon;

    g = 9810;
    ankle_marker = [side, 'ANK'];
    lz = abs(markers.(ankle_marker)(1:end-1, 3) - COM_Z(1:end-1));
    k  = sqrt(g ./ (lz + 1e-6));
    xCOM_AP = COM_AP(1:end-1) + vAP ./ k;
    xCOM_ML = COM_ML(1:end-1) + vML ./ k;

    heel_marker = [side, 'HEE'];
    deltaY = markers.(heel_marker)(to, 2) - markers.(heel_marker)(hs, 2);
    if deltaY > 0, dirAP = 1; else, dirAP = -1; end

    m5_side = [side, 'M5']; m5_opp = [oppside, 'M5'];
    if (markers.(m5_side)(1,1) - markers.(m5_opp)(1,1)) < 0, dirML = -1; else, dirML = 1; end

    toe_marker = [side, 'TOE'];
    HEEL_AP  = markers.(heel_marker)(:, 2);
    TOE_AP   = markers.(toe_marker)(:, 2);
    ANKLE_ML = markers.(ankle_marker)(:, 1);
    M5_ML    = markers.(m5_side)(:, 1);

    mos_ap_heel = dirAP * (HEEL_AP(hs:ho) - xCOM_AP(hs:ho));
    mos_ap_toe  = dirAP * (TOE_AP(ho:to)  - xCOM_AP(ho:to));
    mos_ap_full = [mos_ap_heel; mos_ap_toe];

    mos_ml_ankle = (ANKLE_ML(hs:ho) - xCOM_ML(hs:ho)) * dirML;
    mos_ml_m5    = (M5_ML(ho:to)   - xCOM_ML(ho:to))  * dirML;
    mos_ml_full  = [mos_ml_ankle; mos_ml_m5];

    x_original = linspace(0, 1, length(mos_ap_full));
    x_norm = linspace(0, 1, n_points);

    curve.mos_ap = interp1(x_original, mos_ap_full, x_norm);
    curve.mos_ml = interp1(x_original, mos_ml_full, x_norm);

    if exist('L0_mm','var') && isnumeric(L0_mm) && isfinite(L0_mm) && L0_mm > 0
        curve.mos_ap_P = 100 * interp1(x_original, mos_ap_full./L0_mm, x_norm);
        curve.mos_ml_P = 100 * interp1(x_original, mos_ml_full./L0_mm, x_norm);
    else
        curve.mos_ap_P = nan(1, n_points);
        curve.mos_ml_P = nan(1, n_points);
        warning('L0_mm invalide ou manquant: les séries %%L0 sont NaN.');
    end

    btkCloseAcquisition(data);
end