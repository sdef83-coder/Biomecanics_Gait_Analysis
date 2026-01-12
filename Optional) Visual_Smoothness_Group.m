%  VISUALISATION SMOOTHNESS (SPARC/LDLJ) PAR GROUPE D'ÂGE
%  - Agrégation des fichiers Smoothness_TrialBased_XXX.mat
%  - 3 partie de Figures séparées : (i) (COM vs STERN) × (SPARC vs LDLJ);
%  (ii) Comparaison COM vs STERN; (iii) Comparaison LDLJ et SPARC
%  - Barplots avec 4 groupes d'âge par condition (Plat/Medium/High)
%  - Sorties : 8 figures au total (i = 4, ii = 2, iii = 2)

clear; clc; close all;

%% === PATHS GÉNÉRAUX ===
root_res_all = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL';
cd(root_res_all);

addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script'));

% Dossier où sont les Smoothness_TrialBased_XX.mat
smooth_res_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Smoothness';

% Dossier de sauvegarde des figures
save_fig_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Smoothness\Fig';
if ~exist(save_fig_dir, 'dir'); mkdir(save_fig_dir); end

%% === GROUPES DE PARTICIPANTS ===
ParticipantGroup;
group_names = fieldnames(Group);

%% === PARAMÈTRES ===
surfaces = {'Plat', 'Medium', 'High'};

% Couleurs par groupe d'âge (même palette que MoS)
color_groups = [
    0    0.6   1;    % JeunesEnfants - Bleu clair
    1    0.5   0;    % Enfants - Orange
    0    0.7   0.3;  % Adolescents - Vert
    0.5  0     1];   % Adultes - Violet

% Métriques à analyser
directions = {'AP', 'ML', 'V', 'Magnitude'};
indices = {'SPARC', 'LDLJ'};
segments = {'COM', 'STERN'};

%% === STRUCTURE POUR STOCKER LES DONNÉES ===
% smoothData.(Segment).(Indice).(Direction).(Groupe).(Surface) = [valeurs]
smoothData = struct();

fprintf('🔄 DÉBUT DE LA COLLECTE DES DONNÉES SMOOTHNESS...\n');

%% === COLLECTE DES DONNÉES PAR GROUPE ===
for g = 1:numel(group_names)
    gname = group_names{g};
    sujets = Group.(gname);
    fprintf('\n===== GROUPE : %s =====\n', gname);
    
    for iS = 1:numel(sujets)
        sujet_id = sujets{iS};
        fprintf('  → Sujet %s\n', sujet_id);
        
        % Charger le fichier Smoothness
        mat_file = fullfile(smooth_res_dir, sprintf('Smoothness_TrialBased_%s.mat', sujet_id));
        
        if ~isfile(mat_file)
            fprintf('    ⚠️ Fichier Smoothness introuvable pour %s, ignoré\n', sujet_id);
            continue;
        end
        
        try
            tmp = load(mat_file, 'results');
            results = tmp.results;
        catch ME
            fprintf('    ⚠️ Erreur chargement %s: %s\n', sujet_id, ME.message);
            continue;
        end
        
        % Extraction des données par surface
        for s = 1:numel(surfaces)
            surf = surfaces{s};
            idx_surf = strcmp(results.Surface, surf);
            res_surf = results(idx_surf, :);
            
            if isempty(res_surf)
                continue;
            end
            
            % Pour chaque segment (COM/STERN)
            for seg_idx = 1:numel(segments)
                segment = segments{seg_idx};
                
                % Pour chaque indice (SPARC/LDLJ)
                for ind_idx = 1:numel(indices)
                    indice = indices{ind_idx};
                    
                    % Pour chaque direction
                    for dir_idx = 1:numel(directions)
                        direction = directions{dir_idx};
                        
                        % Nom de la colonne dans la table
                        col_name = sprintf('%s_%s_%s', segment, indice, direction);
                        
                        if ~ismember(col_name, res_surf.Properties.VariableNames)
                            continue;
                        end
                        
                        % Extraire les valeurs (moyenne sur les essais du sujet pour cette surface)
                        values = res_surf.(col_name);
                        clean_values = values(~isnan(values));
                        
                        if isempty(clean_values)
                            continue;
                        end
                        
                        % Moyenne par sujet (pour ne pas sur-représenter les sujets avec plus d'essais)
                        subj_mean = mean(clean_values);
                        
                        % Initialiser la structure si nécessaire
                        if ~isfield(smoothData, segment)
                            smoothData.(segment) = struct();
                        end
                        if ~isfield(smoothData.(segment), indice)
                            smoothData.(segment).(indice) = struct();
                        end
                        if ~isfield(smoothData.(segment).(indice), direction)
                            smoothData.(segment).(indice).(direction) = struct();
                        end
                        if ~isfield(smoothData.(segment).(indice).(direction), gname)
                            smoothData.(segment).(indice).(direction).(gname) = struct();
                        end
                        if ~isfield(smoothData.(segment).(indice).(direction).(gname), surf)
                            smoothData.(segment).(indice).(direction).(gname).(surf) = [];
                        end
                        
                        % Ajouter la moyenne du sujet
                        smoothData.(segment).(indice).(direction).(gname).(surf) = ...
                            [smoothData.(segment).(indice).(direction).(gname).(surf); subj_mean];
                    end
                end
            end
        end
    end
end

fprintf('\n✅ COLLECTE TERMINÉE.\n');

%% === CALCUL DES STATISTIQUES DE GROUPE ===
fprintf('\n📊 Calcul des statistiques par groupe...\n');

groupStats = struct();

for seg_idx = 1:numel(segments)
    segment = segments{seg_idx};
    
    for ind_idx = 1:numel(indices)
        indice = indices{ind_idx};
        
        for dir_idx = 1:numel(directions)
            direction = directions{dir_idx};
            
            for g = 1:numel(group_names)
                gname = group_names{g};
                
                for s = 1:numel(surfaces)
                    surf = surfaces{s};
                    
                    % Vérifier si les données existent
                    if isfield(smoothData, segment) && ...
                       isfield(smoothData.(segment), indice) && ...
                       isfield(smoothData.(segment).(indice), direction) && ...
                       isfield(smoothData.(segment).(indice).(direction), gname) && ...
                       isfield(smoothData.(segment).(indice).(direction).(gname), surf)
                        
                        data = smoothData.(segment).(indice).(direction).(gname).(surf);
                        
                        if ~isempty(data)
                            % Initialiser la structure
                            if ~isfield(groupStats, segment)
                                groupStats.(segment) = struct();
                            end
                            if ~isfield(groupStats.(segment), indice)
                                groupStats.(segment).(indice) = struct();
                            end
                            if ~isfield(groupStats.(segment).(indice), direction)
                                groupStats.(segment).(indice).(direction) = struct();
                            end
                            if ~isfield(groupStats.(segment).(indice).(direction), gname)
                                groupStats.(segment).(indice).(direction).(gname) = struct();
                            end
                            
                            % Calculer les statistiques
                            groupStats.(segment).(indice).(direction).(gname).(surf).mean = mean(data);
                            groupStats.(segment).(indice).(direction).(gname).(surf).std = std(data);
                            groupStats.(segment).(indice).(direction).(gname).(surf).n = numel(data);
                            groupStats.(segment).(indice).(direction).(gname).(surf).sem = std(data) / sqrt(numel(data));
                        end
                    end
                end
            end
        end
    end
end

fprintf('✅ Statistiques calculées.\n');

%% === SAUVEGARDE DES DONNÉES ===
save_data_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Smoothness\Groupes';
if ~exist(save_data_dir, 'dir'); mkdir(save_data_dir); end

smoothnessGroupData = struct();
smoothnessGroupData.smoothData = smoothData;
smoothnessGroupData.groupStats = groupStats;
smoothnessGroupData.surfaces = surfaces;
smoothnessGroupData.groupNames = group_names;
smoothnessGroupData.directions = directions;
smoothnessGroupData.indices = indices;
smoothnessGroupData.segments = segments;
smoothnessGroupData.colors = color_groups;

save(fullfile(save_data_dir, 'Smoothness_groupes.mat'), 'smoothnessGroupData');
fprintf('💾 Données sauvegardées dans : %s\\Smoothness_groupes.mat\n', save_data_dir);

%% ================================================================
%  GÉNÉRATION DES FIGURES - BARPLOTS PAR SEGMENT
%  4 FIGURES : (COM vs STERN) × (SPARC vs LDLJ)
%  Chaque figure : 4 subplots (AP, ML, V, Magnitude)
%  Barres = Moyenne de groupe ; Erreurs = SEM ; Points = sujets (plus transparents)
% ================================================================
fprintf('\n🎨 Génération des figures (4 subplots par figure)...\n');

rng(0); % jitter reproductible

% Paramètres points participants
show_points        = true;
point_size         = 30;
point_face_alpha   = 0.35;   % <<< moins opaque (à ajuster si besoin)
point_edge_alpha   = 0.25;   % <<< contour très léger
point_edge_width   = 0.3;
jitter_frac        = 0.35;   % dispersion horizontale des points

% Paramètres barres
bar_width = 0.18;

for seg_idx = 1:numel(segments)
    segment = segments{seg_idx};

    for ind_idx = 1:numel(indices)
        indice = indices{ind_idx};

        % Créer une figure unique pour ce (segment, indice)
        fig = figure('Name', sprintf('%s - %s', segment, indice), ...
                     'Position', [50, 50, 1400, 900], ...
                     'Color', 'w');

        for dir_idx = 1:numel(directions)
            direction = directions{dir_idx};

            subplot(2,2,dir_idx);
            hold on; box on; grid on;

            % Préparer matrices stats (groupStats) pour ce subplot
            means_matrix = nan(numel(group_names), numel(surfaces));
            sems_matrix  = nan(numel(group_names), numel(surfaces));

            for g = 1:numel(group_names)
                gname = group_names{g};
                for s = 1:numel(surfaces)
                    surf = surfaces{s};

                    if isfield(groupStats, segment) && ...
                       isfield(groupStats.(segment), indice) && ...
                       isfield(groupStats.(segment).(indice), direction) && ...
                       isfield(groupStats.(segment).(indice).(direction), gname) && ...
                       isfield(groupStats.(segment).(indice).(direction).(gname), surf)

                        means_matrix(g, s) = groupStats.(segment).(indice).(direction).(gname).(surf).mean;
                        sems_matrix(g, s)  = groupStats.(segment).(indice).(direction).(gname).(surf).sem;
                    end
                end
            end

            % Si subplot vide (tout NaN) : afficher quand même axes mais sans plot
            if all(isnan(means_matrix(:)))
                title(sprintf('%s - %s', indice, direction), 'FontSize', 12, 'FontWeight', 'bold');
                set(gca, 'XTick', 1:numel(surfaces), 'XTickLabel', surfaces);
                xlabel('Surface', 'FontSize', 11, 'FontWeight', 'bold');
                ylabel(sprintf('%s %s (ua)', indice, direction), 'FontSize', 11, 'FontWeight', 'bold');
                continue;
            end

            % Positions
            x = 1:numel(surfaces);

            % Barres + SEM + points
            for g = 1:numel(group_names)
                gname  = group_names{g};
                offset = (g - 2.5) * bar_width;

                % Barres
                bar(x + offset, means_matrix(g, :), bar_width, ...
                    'FaceColor', color_groups(g, :), ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.7, ...
                    'DisplayName', gname);

                % SEM
                for s = 1:numel(surfaces)
                    if ~isnan(means_matrix(g, s)) && ~isnan(sems_matrix(g, s))
                        errorbar(x(s) + offset, means_matrix(g, s), sems_matrix(g, s), ...
                            'k', 'LineStyle', 'none', 'LineWidth', 1.2, ...
                            'CapSize', 5, 'HandleVisibility', 'off');
                    end
                end

                % Points individuels (sujets)
                if show_points
                    for s = 1:numel(surfaces)
                        surf = surfaces{s};

                        if isfield(smoothData, segment) && ...
                           isfield(smoothData.(segment), indice) && ...
                           isfield(smoothData.(segment).(indice), direction) && ...
                           isfield(smoothData.(segment).(indice).(direction), gname) && ...
                           isfield(smoothData.(segment).(indice).(direction).(gname), surf)

                            data = smoothData.(segment).(indice).(direction).(gname).(surf);
                            data = data(~isnan(data));

                            if ~isempty(data)
                                jitter = (rand(size(data)) - 0.5) * (bar_width * 2 * jitter_frac);
                                xpts = (x(s) + offset) + jitter;

                                sc = scatter(xpts, data, point_size, ...
                                    'MarkerFaceColor', color_groups(g, :), ...
                                    'MarkerEdgeColor', 'k', ...
                                    'LineWidth', point_edge_width, ...
                                    'MarkerFaceAlpha', point_face_alpha, ...
                                    'MarkerEdgeAlpha', point_edge_alpha);

                                sc.HandleVisibility = 'off';
                            end
                        end
                    end
                end
            end

            % Axes / labels
            set(gca, 'XTick', x, 'XTickLabel', surfaces);
            xlabel('Surface', 'FontSize', 11, 'FontWeight', 'bold');

            if strcmp(indice, 'SPARC')
                ylabel_text = sprintf('%s SPARC %s (ua)', segment, direction);
            else
                ylabel_text = sprintf('%s LDLJ %s (ua)', segment, direction);
            end
            ylabel(ylabel_text, 'FontSize', 11, 'FontWeight', 'bold');

            title(sprintf('%s - %s', indice, direction), 'FontSize', 12, 'FontWeight', 'bold');

            % Limites Y
            yl = ylim;
            ylim([yl(1), yl(2) * 1.15]);

            % Légende seulement sur le 1er subplot
            if dir_idx == 1
                legend('Location', 'best', 'FontSize', 9);
            end

            hold off;
        end

        sgtitle(sprintf('%s - %s (Moyenne ± SEM + sujets)', segment, indice), ...
                'FontSize', 16, 'FontWeight', 'bold');

        % Sauvegarde
        filename = sprintf('BARPLOT_%s_%s_4subplots.png', segment, indice);
        saveas(fig, fullfile(save_fig_dir, filename));
        fprintf('  ✅ %s\n', filename);

        close(fig);
    end
end

%% ================================================================
%  FIGURES COMPARATIVES : COM vs STERN (CÔTE À CÔTE) - CORRIGÉ (anti-chevauchement)
% ================================================================
fprintf('\n🎨 Génération des figures comparatives COM vs STERN...\n');

for ind_idx = 1:numel(indices)
    indice = indices{ind_idx};

    % Figure avec 4 subplots (une par direction)
    fig = figure('Name', sprintf('Comparaison COM vs STERN - %s', indice), ...
                 'Position', [50, 50, 1600, 1000], ...
                 'Color', 'w');

    for dir_idx = 1:numel(directions)
        direction = directions{dir_idx};

        subplot(2, 2, dir_idx);
        hold on; box on; grid on;

        % Positions des surfaces
        x = 1:numel(surfaces);

        % --------- PARAMÈTRES ANTI-CHEVAUCHEMENT ----------
        nGroups   = numel(group_names);
        nSegments = 2;                    % COM + STERN
        nBars     = nGroups * nSegments;  % 8 barres par surface

        total_width = 0.80;               % largeur totale occupée autour de chaque surface (<= 1 recommandé)
        bar_width   = total_width / nBars;

        % Bord gauche du "pack" de barres pour chaque surface
        x0 = x - total_width/2;
        % --------------------------------------------------

        % Tracer COM et STERN pour chaque groupe
        for g = 1:nGroups
            gname = group_names{g};

            % Préparer les données COM
            means_com = nan(1, numel(surfaces));
            sems_com  = nan(1, numel(surfaces));

            for s = 1:numel(surfaces)
                surf = surfaces{s};
                if isfield(groupStats, 'COM') && ...
                   isfield(groupStats.COM, indice) && ...
                   isfield(groupStats.COM.(indice), direction) && ...
                   isfield(groupStats.COM.(indice).(direction), gname) && ...
                   isfield(groupStats.COM.(indice).(direction).(gname), surf)

                    means_com(s) = groupStats.COM.(indice).(direction).(gname).(surf).mean;
                    sems_com(s)  = groupStats.COM.(indice).(direction).(gname).(surf).sem;
                end
            end

            % Préparer les données STERN
            means_stern = nan(1, numel(surfaces));
            sems_stern  = nan(1, numel(surfaces));

            for s = 1:numel(surfaces)
                surf = surfaces{s};
                if isfield(groupStats, 'STERN') && ...
                   isfield(groupStats.STERN, indice) && ...
                   isfield(groupStats.STERN.(indice), direction) && ...
                   isfield(groupStats.STERN.(indice).(direction), gname) && ...
                   isfield(groupStats.STERN.(indice).(direction).(gname), surf)

                    means_stern(s) = groupStats.STERN.(indice).(direction).(gname).(surf).mean;
                    sems_stern(s)  = groupStats.STERN.(indice).(direction).(gname).(surf).sem;
                end
            end

            % Indices de barres dans le pack (1..nBars) : COM puis STERN pour chaque groupe
            k_com   = (g-1)*nSegments + 1;   % 1,3,5,7
            k_stern = (g-1)*nSegments + 2;   % 2,4,6,8

            % Positions x exactes (centres des barres)
            x_com   = x0 + (k_com   - 0.5) * bar_width;
            x_stern = x0 + (k_stern - 0.5) * bar_width;

            % Barres COM (plus sombres)
            bar(x_com, means_com, bar_width, ...
                'FaceColor', min(color_groups(g, :) * 0.7, 1), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.8, ...
                'DisplayName', sprintf('%s COM', gname));

            % Barres STERN (plus claires)
            bar(x_stern, means_stern, bar_width, ...
                'FaceColor', color_groups(g, :), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.5, ...
                'DisplayName', sprintf('%s STERN', gname));

            % Barres d'erreur COM
            for s = 1:numel(surfaces)
                if ~isnan(means_com(s)) && ~isnan(sems_com(s))
                    errorbar(x_com(s), means_com(s), sems_com(s), ...
                             'k', 'LineStyle', 'none', 'LineWidth', 1.2, ...
                             'CapSize', 4, 'HandleVisibility', 'off');
                end
            end

            % Barres d'erreur STERN
            for s = 1:numel(surfaces)
                if ~isnan(means_stern(s)) && ~isnan(sems_stern(s))
                    errorbar(x_stern(s), means_stern(s), sems_stern(s), ...
                             'k', 'LineStyle', 'none', 'LineWidth', 1.2, ...
                             'CapSize', 4, 'HandleVisibility', 'off');
                end
            end
        end

        % Personnalisation
        set(gca, 'XTick', x, 'XTickLabel', surfaces);
        xlim([0.5, numel(surfaces)+0.5]); % évite les débordements visuels
        xlabel('Surface', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel(sprintf('%s %s (ua)', indice, direction), 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s - %s', indice, direction), 'FontSize', 12, 'FontWeight', 'bold');

        % Légende seulement sur le premier subplot
        if dir_idx == 1
            legend('Location', 'eastoutside', 'FontSize', 8);
        end

        ylims = ylim;
        ylim([ylims(1), ylims(2) * 1.15]);

        hold off;
    end

    sgtitle(sprintf('Comparaison COM vs STERNUM - %s (Moyenne ± SEM)', indice), ...
            'FontSize', 16, 'FontWeight', 'bold');

    % Sauvegarder
    filename = sprintf('Comparison_COM_vs_STERN_%s.png', indice);
    saveas(fig, fullfile(save_fig_dir, filename));
    fprintf('  ✅ %s\n', filename);

    close(fig);
end

%% ================================================================
%  FIGURES COMPARATIVES : LDLJ vs SPARC (CÔTE À CÔTE) - ANTI-CHEVAUCHEMENT
%  - Une figure par segment (COM / STERN)
%  - 4 subplots (AP, ML, V, Magnitude)
%  - Barres : 4 groupes × 2 indices (SPARC, LDLJ) = 8 barres par surface
% ================================================================
fprintf('\n🎨 Génération des figures comparatives LDLJ vs SPARC...\n');

for seg_idx = 1:numel(segments)
    segment = segments{seg_idx};

    % Figure avec 4 subplots (une par direction)
    fig = figure('Name', sprintf('Comparaison LDLJ vs SPARC - %s', segment), ...
                 'Position', [50, 50, 1600, 1000], ...
                 'Color', 'w');

    for dir_idx = 1:numel(directions)
        direction = directions{dir_idx};

        subplot(2, 2, dir_idx);
        hold on; box on; grid on;

        % Positions des surfaces
        x = 1:numel(surfaces);

        % --------- PARAMÈTRES ANTI-CHEVAUCHEMENT ----------
        nGroups   = numel(group_names);
        nIndices  = 2;                    % SPARC + LDLJ
        nBars     = nGroups * nIndices;   % 8 barres par surface

        total_width = 0.80;               % largeur totale occupée autour de chaque surface (<= 1 recommandé)
        bar_width   = total_width / nBars;

        % Bord gauche du "pack" de barres pour chaque surface
        x0 = x - total_width/2;
        % --------------------------------------------------

        % Tracer SPARC et LDLJ pour chaque groupe
        for g = 1:nGroups
            gname = group_names{g};

            % Préparer les données SPARC
            means_sparc = nan(1, numel(surfaces));
            sems_sparc  = nan(1, numel(surfaces));

            for s = 1:numel(surfaces)
                surf = surfaces{s};
                if isfield(groupStats, segment) && ...
                   isfield(groupStats.(segment), 'SPARC') && ...
                   isfield(groupStats.(segment).SPARC, direction) && ...
                   isfield(groupStats.(segment).SPARC.(direction), gname) && ...
                   isfield(groupStats.(segment).SPARC.(direction).(gname), surf)

                    means_sparc(s) = groupStats.(segment).SPARC.(direction).(gname).(surf).mean;
                    sems_sparc(s)  = groupStats.(segment).SPARC.(direction).(gname).(surf).sem;
                end
            end

            % Préparer les données LDLJ
            means_ldlj = nan(1, numel(surfaces));
            sems_ldlj  = nan(1, numel(surfaces));

            for s = 1:numel(surfaces)
                surf = surfaces{s};
                if isfield(groupStats, segment) && ...
                   isfield(groupStats.(segment), 'LDLJ') && ...
                   isfield(groupStats.(segment).LDLJ, direction) && ...
                   isfield(groupStats.(segment).LDLJ.(direction), gname) && ...
                   isfield(groupStats.(segment).LDLJ.(direction).(gname), surf)

                    means_ldlj(s) = groupStats.(segment).LDLJ.(direction).(gname).(surf).mean;
                    sems_ldlj(s)  = groupStats.(segment).LDLJ.(direction).(gname).(surf).sem;
                end
            end

            % Indices de barres dans le pack (1..nBars) : SPARC puis LDLJ pour chaque groupe
            k_sparc = (g-1)*nIndices + 1;   % 1,3,5,7
            k_ldlj  = (g-1)*nIndices + 2;   % 2,4,6,8

            % Positions x exactes (centres des barres)
            x_sparc = x0 + (k_sparc - 0.5) * bar_width;
            x_ldlj  = x0 + (k_ldlj  - 0.5) * bar_width;

            % Barres SPARC (plus sombres)
            bar(x_sparc, means_sparc, bar_width, ...
                'FaceColor', min(color_groups(g, :) * 0.7, 1), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.8, ...
                'DisplayName', sprintf('%s SPARC', gname));

            % Barres LDLJ (plus claires)
            bar(x_ldlj, means_ldlj, bar_width, ...
                'FaceColor', color_groups(g, :), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.5, ...
                'DisplayName', sprintf('%s LDLJ', gname));

            % Erreurs SPARC
            for s = 1:numel(surfaces)
                if ~isnan(means_sparc(s)) && ~isnan(sems_sparc(s))
                    errorbar(x_sparc(s), means_sparc(s), sems_sparc(s), ...
                             'k', 'LineStyle', 'none', 'LineWidth', 1.2, ...
                             'CapSize', 4, 'HandleVisibility', 'off');
                end
            end

            % Erreurs LDLJ
            for s = 1:numel(surfaces)
                if ~isnan(means_ldlj(s)) && ~isnan(sems_ldlj(s))
                    errorbar(x_ldlj(s), means_ldlj(s), sems_ldlj(s), ...
                             'k', 'LineStyle', 'none', 'LineWidth', 1.2, ...
                             'CapSize', 4, 'HandleVisibility', 'off');
                end
            end
        end

        % Personnalisation
        set(gca, 'XTick', x, 'XTickLabel', surfaces);
        xlim([0.5, numel(surfaces)+0.5]);
        xlabel('Surface', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel(sprintf('%s %s (ua)', segment, direction), 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s - %s', segment, direction), 'FontSize', 12, 'FontWeight', 'bold');

        % Légende seulement sur le premier subplot
        if dir_idx == 1
            legend('Location', 'eastoutside', 'FontSize', 8);
        end

        ylims = ylim;
        ylim([ylims(1), ylims(2) * 1.15]);

        hold off;
    end

    sgtitle(sprintf('Comparaison LDLJ vs SPARC - %s (Moyenne ± SEM)', segment), ...
            'FontSize', 16, 'FontWeight', 'bold');

    % Sauvegarder
    filename = sprintf('Comparison_LDLJ_vs_SPARC_%s.png', segment);
    saveas(fig, fullfile(save_fig_dir, filename));
    fprintf('  ✅ %s\n', filename);

    close(fig);
end

fprintf('\n🎯 TERMINÉ! Toutes les figures sont dans : %s\n', save_fig_dir);