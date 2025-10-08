%% CLUSTERING NON SUPERVISE (K-MEANS) SUR LES VARIABLES SPATIO-TEMPORELLES A LA MARCHE
clc; clear; close all;

% === Dossier d'E/S ===
root_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result';
cd(root_path)
load('SpatioTemporalDATA.mat')

outdir = fullfile(root_path, 'Fig', 'Clustering');
if ~exist(outdir, 'dir'); mkdir(outdir); end
ts = string(datetime('now','Format','yyyyMMdd_HHmm')); % horodatage

% 1) Variables d'intérêt (Mean + CV)
varsMean = { ...
    'Mean_Double support time (%)', ...
    'Mean_Stride width (mm)', ...
    'Mean_Gait speed (m.s^{-1})', ...
    'Mean_Normalized Walk ratio (ua)', ...
    'Mean_Normalized Step length (ua)', ...
    'Mean_Normalized Cadence (ua)'};
varsCV = strrep(varsMean, 'Mean_', 'CV_');
allVars = [varsMean, varsCV];

% 2) Conditions et couleurs associées
conds = {'Plat','Medium','High'};
colors = struct('Plat','b','Medium','g','High','r'); %#ok<NASGU>

% 3) Préallocation
dataMat = [];
meta = table();  % Condition, Participant, AgeMonths

% 4) Lecture et concaténation des données
for iC = 1:numel(conds)
    cond = conds{iC};
    T = SpatioTemporalDATA.ALL.(cond);
    mask = all(~ismissing(T(:, allVars)), 2);
    Tsel = T(mask, :);
    X = Tsel{:, allVars};
    dataMat = [dataMat; X];
    n = size(X,1);
    meta = [meta; table(repmat(string(cond),n,1), Tsel.Participant, Tsel.AgeMonths, ...
        'VariableNames',{'Condition','Participant','AgeMonths'})];
end

% 4b) Ajoute colonne Groupe d'âge (catégories fixes)
meta.AgeGroup = arrayfun(@(m) age_group_from_months(m), meta.AgeMonths, 'UniformOutput', false);
meta.AgeGroup = string(meta.AgeGroup);

% 5) Normalisation (z-score)
dataNorm = zscore(dataMat);

% 6) ACP (sur les variables normalisées) — BASE COMMUNE
[coeff, score, ~, ~, explained] = pca(dataNorm);

% === (A) ANALYSE GLOBALE (toutes conditions) ============================

% --- Biplot PC1-PC2 + cercle de corrélation
fig_biplot = figure('Color','w');
biplot(coeff(:,1:2), 'Scores', score(:,1:2), 'VarLabels', allVars);
axis equal; hold on;
theta = linspace(0,2*pi,100); plot(cos(theta), sin(theta), 'k--','LineWidth',1);
title('ACP : PC1 & PC2 avec cercle de corrélation');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
grid on;
exportgraphics(gca, fullfile(outdir, "01_Biplot_PC1_PC2_"+ts+".png"), 'Resolution',300);

% --- Scatter PC1-PC2 : surfaces & tranches d'âge
fig_scatter = figure('Color','w'); hold on;
ageEdges = [0,72,144,216,inf];  % mois: 0-6, 6-12, 12-18, >18 ans
ageLabelsAge = {'0-6 ans','6-12 ans','12-18 ans','>18 ans'};
for iC = 1:numel(conds)
    cond = conds{iC}; 
    idxC = meta.Condition==cond;
    for iA = 1:numel(ageLabelsAge)
        idxA = meta.AgeMonths>ageEdges(iA) & meta.AgeMonths<=ageEdges(iA+1);
        idx = idxC & idxA;
        if any(idx)
            scatter(score(idx,1), score(idx,2), 50, 'filled', ...
                'DisplayName', sprintf('%s | %s', cond, ageLabelsAge{iA}));
        end
    end
end
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('Observations : Surfaces & Tranches d age');
legend('Location','bestoutside'); grid on;
exportgraphics(gca, fullfile(outdir, "02_Scatter_PC1_PC2_AgeSurface_"+ts+".png"), 'Resolution',300);

% --- Choix automatique de k : méthode du coude (global)
kmax = 10; rng(42);
sumd_global = elbow_sumd(score(:,1:2), kmax);
[k_global, ~] = elbow_k(sumd_global);

fig_elbow = figure('Color','w');
plot(1:kmax, sumd_global, '-o', 'LineWidth',1.5); hold on;
plot([1 kmax], [sumd_global(1) sumd_global(end)], 'k--');
plot(k_global, sumd_global(k_global), 'rp', 'MarkerFaceColor','r', 'MarkerSize',12);
xlabel('Nombre de clusters k'); ylabel('Somme des distances intra-cluster');
title(sprintf('Méthode du coude — k* = %d (global)', k_global));
grid on;
exportgraphics(gca, fullfile(outdir, "03_Elbow_GLOBAL_"+ts+".png"), 'Resolution',300);

% --- Clustering global
[idxCluster_global, Cc_global] = kmeans(score(:,1:2), k_global, 'Replicates',20);

% --- Visualisation clusters (global) + ANNOTATIONS ÂGE
fig_clusters = figure('Color','w'); hold on;
colorsK = lines(k_global);
for i = 1:k_global
    scatter(score(idxCluster_global==i,1), score(idxCluster_global==i,2), 60, colorsK(i,:), 'filled', ...
        'DisplayName', sprintf('Cluster %d', i));
end
plot(Cc_global(:,1), Cc_global(:,2), 'kx', 'MarkerSize',15, 'LineWidth',2, 'DisplayName','Centres');
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title(sprintf('Clustering k-means (k=%d) sur PC1-PC2 — GLOBAL', k_global));
grid on; legend('Location','bestoutside');

% >>> Annotations âge moyen ± SD sur la figure globale
AgeStats_global = annotate_age_on_clusters(gca, score(:,1:2), idxCluster_global, Cc_global, meta.AgeMonths);

% --- Sauvegardes liées au global
exportgraphics(gca, fullfile(outdir, "04_Clusters_GLOBAL_k"+k_global+"_"+ts+".png"), 'Resolution',300);

% Profils moyens (z-score) par cluster (global)
meanVars_global = array2table(zeros(k_global, numel(allVars)), 'VariableNames', allVars, ...
    'RowNames', arrayfun(@(i)sprintf('Cluster_%d',i),1:k_global,'UniformOutput',false));
for i = 1:k_global
    meanVars_global{i, :} = mean(dataNorm(idxCluster_global==i, :), 1, 'omitnan');
end

% === HEATMAPS ANALYSE GLOBALE ===

% --- Heatmap globale : profils moyens z-score par cluster ---
fig_heatmap_global = figure('Color','w');
fig_heatmap_global.Position = [100, 100, 1200, 600];

% Matrice des profils moyens (clusters x variables)
heatmap_data_global = table2array(meanVars_global);

% Créer la heatmap
h_global = heatmap(allVars, meanVars_global.Properties.RowNames, heatmap_data_global, ...
    'Title', sprintf('Profils moyens z-score par cluster (Global, k=%d)', k_global), ...
    'XLabel', 'Variables spatio-temporelles', ...
    'YLabel', 'Clusters', ...
    'Colormap', blueyellow(256), ...
    'ColorLimits', [-max(abs(heatmap_data_global(:))), max(abs(heatmap_data_global(:)))]);

% Personnaliser l'affichage
h_global.FontSize = 10;
h_global.CellLabelFormat = '%.2f';
h_global.GridVisible = 'on';

% Augmenter la taille du titre
try
    h_global.Title.FontSize = 16;
    h_global.Title.FontWeight = 'bold';
catch
    warning('Impossible de modifier la taille du titre');
end

% Rotation des étiquettes X pour meilleure lisibilité
try
    % Pour les versions récentes de MATLAB
    h_global.XDisplayLabels = cellfun(@(x) strrep(x, ' ', {' '; ''}), allVars, 'UniformOutput', false);
catch
    % Méthode alternative si la propriété n'existe pas
    warning('Impossible de modifier les étiquettes X de la heatmap');
end

% Sauvegarder
exportgraphics(h_global, fullfile(outdir, "05_Heatmap_Profils_GLOBAL_k"+k_global+"_"+ts+".png"), 'Resolution',300);

% --- Heatmap globale avec annotations statistiques ---
fig_heatmap_stats_global = figure('Color','w');
fig_heatmap_stats_global.Position = [100, 150, 1400, 700];

% Créer une version annotée avec informations sur les clusters
cluster_info_global = cell(k_global, 1);
for i = 1:k_global
    n_cluster = sum(idxCluster_global == i);
    age_mean = AgeStats_global.AgeMonths_mean(i);
    cluster_info_global{i} = sprintf('Cluster %d (n=%d, %.1f ans)', i, n_cluster, age_mean/12);
end

h_stats_global = heatmap(allVars, cluster_info_global, heatmap_data_global, ...
    'Title', sprintf('Profils moyens z-score avec informations clusters (Global, k=%d)', k_global), ...
    'XLabel', 'Variables spatio-temporelles', ...
    'YLabel', 'Clusters (taille, âge moyen)', ...
    'Colormap', blueyellow(256), ...
    'ColorLimits', [-max(abs(heatmap_data_global(:))), max(abs(heatmap_data_global(:)))]);

h_stats_global.FontSize = 9;
h_stats_global.CellLabelFormat = '%.2f';
h_stats_global.GridVisible = 'on';

% Augmenter la taille du titre
try
    h_stats_global.Title.FontSize = 16;
    h_stats_global.Title.FontWeight = 'bold';
catch
    warning('Impossible de modifier la taille du titre');
end

% Rotation des étiquettes X
try
    % Pour les versions récentes de MATLAB
    h_stats_global.XDisplayLabels = cellfun(@(x) strrep(x, ' ', {' '; ''}), allVars, 'UniformOutput', false);
catch
    % Méthode alternative si la propriété n'existe pas
    warning('Impossible de modifier les étiquettes X de la heatmap');
end

exportgraphics(h_stats_global, fullfile(outdir, "06_Heatmap_Profils_Stats_GLOBAL_k"+k_global+"_"+ts+".png"), 'Resolution',300);

% CSV global
writetable(AgeStats_global, fullfile(outdir, "Age_Stats_GLOBAL_"+ts+".csv"));
writetable( addvars(meanVars_global,(1:k_global)','Before',1,'NewVariableNames','ClusterID'), ...
            fullfile(outdir, "ProfilMoyen_zscore_GLOBAL_"+ts+".csv"), ...
            'WriteRowNames', true, 'Delimiter',',');
assignTable_global = table(meta.Participant, meta.Condition, meta.AgeMonths, meta.AgeGroup, ...
                           score(:,1), score(:,2), idxCluster_global, ...
    'VariableNames', {'Participant','Condition','AgeMonths','AgeGroup','PC1','PC2','Cluster'});
writetable(assignTable_global, fullfile(outdir, "Assignations_PC12_GLOBAL_"+ts+".csv"));
writetable(table((1:kmax)', sumd_global', 'VariableNames', {'k','SumIntraClusterDistances'}), ...
           fullfile(outdir, "Elbow_GLOBAL_k_SumD_"+ts+".csv"));
writetable(table((1:numel(explained))', explained, 'VariableNames', {'PC','ExplainedVariancePct'}), ...
           fullfile(outdir, "ACP_ExplainedVariance_GLOBAL_"+ts+".csv"));

% ===>>> AFFICHAGE CONSOLE — COMPOSITION DES CLUSTERS (GLOBAL) ------------
fprintf('\n=== COMPOSITION DES CLUSTERS (GLOBAL) ===\n');
print_cluster_composition(meta, idxCluster_global, conds);

% === (B) ANALYSE PAR CONDITION avec HEATMAPS ===
independentPCA = false;  % true = ACP recalculée par condition (pas comparables entre conditions)

for iC = 1:numel(conds)
    cond = conds{iC};
    idxCond = meta.Condition==cond;
    nCond = sum(idxCond);
    if nCond < 3
        warning('Condition %s: trop peu d''observations (%d). Clustering sauté.', cond, nCond);
        continue;
    end

    if independentPCA
        % ACP propre à la condition
        Xnorm_c = dataNorm(idxCond, :);
        [coeff_c, score_c, ~, ~, explained_c] = pca(Xnorm_c); %#ok<ASGLU>
        pcs   = score_c(:,1:2);
        expl1 = explained_c(1); expl2 = explained_c(2);
    else
        % Réutilise l'ACP globale
        pcs   = score(idxCond,1:2);
        expl1 = explained(1); expl2 = explained(2);
        Xnorm_c = dataNorm(idxCond, :); % pour profils moyens
    end

    % Coude conditionnel (k auto)
    kmax_c = min(10, max(2, nCond-1));
    sumd_c = elbow_sumd(pcs, kmax_c);
    [k_c, ~] = elbow_k(sumd_c);

    % K-means sur la condition
    [idxC, CcC] = kmeans(pcs, k_c, 'Replicates',20);

    % Figures — clusters (condition) + ANNOTATIONS ÂGE
    figc = figure('Color','w'); hold on;
    colsC = lines(k_c);
    for ic = 1:k_c
        scatter(pcs(idxC==ic,1), pcs(idxC==ic,2), 60, colsC(ic,:), 'filled', ...
            'DisplayName', sprintf('Cluster %d', ic));
    end
    plot(CcC(:,1), CcC(:,2), 'kx', 'MarkerSize',15, 'LineWidth',2, 'DisplayName','Centres');
    xlabel(sprintf('PC1 (%.1f%%)', expl1));
    ylabel(sprintf('PC2 (%.1f%%)', expl2));
    title(sprintf('Clustering k-means — %s (k=%d)', cond, k_c));
    grid on; legend('Location','bestoutside');

    % >>> Annotations âge moyen ± SD sur la figure conditionnelle
    meta_c = meta(idxCond,:);
    AgeStats_c = annotate_age_on_clusters(gca, pcs, idxC, CcC, meta_c.AgeMonths);

    % Export figure clusters
    exportgraphics(gca, fullfile(outdir, "B1_Clusters_"+cond+"_k"+k_c+"_"+ts+".png"), 'Resolution',300);

    % Figure — coude (condition)
    fig_elb_c = figure('Color','w');
    plot(1:kmax_c, sumd_c, '-o','LineWidth',1.5); hold on;
    plot([1 kmax_c], [sumd_c(1) sumd_c(end)], 'k--');
    plot(k_c, sumd_c(k_c), 'rp', 'MarkerFaceColor','r', 'MarkerSize',12);
    xlabel('k'); ylabel('Somme des distances intra-cluster');
    title(sprintf('Méthode du coude — %s (k*=%d)', cond, k_c));
    grid on;
    exportgraphics(gca, fullfile(outdir, "B0_Elbow_"+cond+"_"+ts+".png"), 'Resolution',300);

    % Profils moyens (z-score) par cluster (condition)
    meanVars_c = array2table(zeros(k_c, numel(allVars)), 'VariableNames', allVars, ...
        'RowNames', arrayfun(@(i)sprintf('Cluster_%d',i),1:k_c,'UniformOutput',false));
    for ic = 1:k_c
        meanVars_c{ic, :} = mean(Xnorm_c(idxC==ic, :), 1, 'omitnan');
    end

    % === HEATMAPS POUR CETTE CONDITION ===
    
    % --- Heatmap simple : profils moyens z-score par cluster ---
    fig_heatmap_cond = figure('Color','w');
    fig_heatmap_cond.Position = [100, 100, 1200, 400 + k_c*50];
    
    heatmap_data_c = table2array(meanVars_c);
    
    h_cond = heatmap(allVars, meanVars_c.Properties.RowNames, heatmap_data_c, ...
        'Title', sprintf('Profils moyens z-score par cluster — %s (k=%d)', cond, k_c), ...
        'XLabel', 'Variables spatio-temporelles', ...
        'YLabel', 'Clusters', ...
        'Colormap', blueyellow(256), ...
        'ColorLimits', [-max(abs(heatmap_data_c(:))), max(abs(heatmap_data_c(:)))]);
    
    h_cond.FontSize = 10;
    h_cond.CellLabelFormat = '%.2f';
    h_cond.GridVisible = 'on';
    
    % Augmenter la taille du titre
    try
        h_cond.Title.FontSize = 16;
        h_cond.Title.FontWeight = 'bold';
    catch
        warning('Impossible de modifier la taille du titre pour %s', cond);
    end
    
    % Rotation des étiquettes X
    try
        % Pour les versions récentes de MATLAB
        h_cond.XDisplayLabels = cellfun(@(x) strrep(x, ' ', {' '; ''}), allVars, 'UniformOutput', false);
    catch
        % Méthode alternative si la propriété n'existe pas
        warning('Impossible de modifier les étiquettes X de la heatmap pour %s', cond);
    end
    
    exportgraphics(h_cond, fullfile(outdir, "B2_Heatmap_Profils_"+cond+"_k"+k_c+"_"+ts+".png"), 'Resolution',300);
    
    % --- Heatmap avec annotations statistiques ---
    fig_heatmap_stats_cond = figure('Color','w');
    fig_heatmap_stats_cond.Position = [100, 150, 1400, 450 + k_c*50];
    
    cluster_info_c = cell(k_c, 1);
    for ic = 1:k_c
        n_cluster_c = sum(idxC == ic);
        age_mean_c = AgeStats_c.AgeMonths_mean(ic);
        cluster_info_c{ic} = sprintf('Cluster %d (n=%d, %.1f ans)', ic, n_cluster_c, age_mean_c/12);
    end
    
    h_stats_cond = heatmap(allVars, cluster_info_c, heatmap_data_c, ...
        'Title', sprintf('Profils moyens z-score avec infos clusters — %s (k=%d)', cond, k_c), ...
        'XLabel', 'Variables spatio-temporelles', ...
        'YLabel', 'Clusters (taille, âge moyen)', ...
        'Colormap', blueyellow(256), ...
        'ColorLimits', [-max(abs(heatmap_data_c(:))), max(abs(heatmap_data_c(:)))]);
    
    h_stats_cond.FontSize = 9;
    h_stats_cond.CellLabelFormat = '%.2f';
    h_stats_cond.GridVisible = 'on';
    
    % Augmenter la taille du titre
    try
        h_stats_cond.Title.FontSize = 16;
        h_stats_cond.Title.FontWeight = 'bold';
    catch
        warning('Impossible de modifier la taille du titre pour %s', cond);
    end
    
    % Rotation des étiquettes X
    try
        % Pour les versions récentes de MATLAB
        h_stats_cond.XDisplayLabels = cellfun(@(x) strrep(x, ' ', {' '; ''}), allVars, 'UniformOutput', false);
    catch
        % Méthode alternative si la propriété n'existe pas
        warning('Impossible de modifier les étiquettes X de la heatmap pour %s', cond);
    end
    
    exportgraphics(h_stats_cond, fullfile(outdir, "B2b_Heatmap_Profils_Stats_"+cond+"_k"+k_c+"_"+ts+".png"), 'Resolution',300);

    % Sauvegardes CSV (condition)
    writetable( addvars(meanVars_c,(1:k_c)','Before',1,'NewVariableNames','ClusterID'), ...
                fullfile(outdir, "B3_ProfilMoyen_zscore_"+cond+"_"+ts+".csv"), ...
                'WriteRowNames', true, 'Delimiter',',');
    assignTable_c = table(meta_c.Participant, meta_c.AgeMonths, meta_c.AgeGroup, ...
                          pcs(:,1), pcs(:,2), idxC, ...
        'VariableNames', {'Participant','AgeMonths','AgeGroup','PC1','PC2','Cluster'});
    writetable(assignTable_c, fullfile(outdir, "B4_Assignations_PC12_"+cond+"_"+ts+".csv"));

    % ===>>> AFFICHAGE CONSOLE — COMPOSITION DES CLUSTERS (PAR CONDITION) ---
    fprintf('\n=== COMPOSITION DES CLUSTERS — Condition %s ===\n', cond);
    % Ici, on affiche uniquement la composition par groupe d'âge (la condition est fixe)
    print_cluster_composition(meta_c, idxC, {cond});  % conds passé pour compatibilité
end

% === HEATMAP COMPARATIVE GLOBALE ===
fprintf('\n=== CRÉATION HEATMAP COMPARATIVE ===\n');

% Concaténer tous les profils moyens
all_profiles = [];
all_labels = {};

% Global d'abord
all_profiles = [all_profiles; table2array(meanVars_global)];
for i = 1:k_global
    all_labels{end+1} = sprintf('Global-C%d (n=%d)', i, sum(idxCluster_global == i));
end

% Puis chaque condition
for iC = 1:numel(conds)
    cond = conds{iC};
    idxCond = meta.Condition==cond;
    nCond = sum(idxCond);
    if nCond < 3, continue; end
    
    % Recalculer les données pour cette condition
    if independentPCA
        Xnorm_c = dataNorm(idxCond, :);
        [~, score_c, ~, ~, ~] = pca(Xnorm_c);
        pcs = score_c(:,1:2);
    else
        pcs = score(idxCond,1:2);
        Xnorm_c = dataNorm(idxCond, :);
    end
    
    kmax_c = min(10, max(2, nCond-1));
    sumd_c = elbow_sumd(pcs, kmax_c);
    [k_c, ~] = elbow_k(sumd_c);
    [idxC, ~] = kmeans(pcs, k_c, 'Replicates',20);
    
    % Profils moyens
    meanVars_c = zeros(k_c, numel(allVars));
    for ic = 1:k_c
        meanVars_c(ic, :) = mean(Xnorm_c(idxC==ic, :), 1, 'omitnan');
        all_labels{end+1} = sprintf('%s-C%d (n=%d)', cond, ic, sum(idxC == ic));
    end
    all_profiles = [all_profiles; meanVars_c];
end

% Créer la heatmap comparative
if ~isempty(all_profiles)
    fig_comparative = figure('Color','w');
    fig_comparative.Position = [100, 50, 1600, 300 + size(all_profiles,1)*30];
    
    h_comp = heatmap(allVars, all_labels, all_profiles, ...
        'Title', 'Comparaison des profils moyens z-score - Toutes conditions', ...
        'XLabel', 'Variables spatio-temporelles', ...
        'YLabel', 'Clusters par condition', ...
        'Colormap', blueyellow(256), ...
        'ColorLimits', [-max(abs(all_profiles(:))), max(abs(all_profiles(:)))]);
    
    h_comp.FontSize = 8;
    h_comp.CellLabelFormat = '%.2f';
    h_comp.GridVisible = 'on';
    
    % Augmenter la taille du titre
    try
        h_comp.Title.FontSize = 16;
        h_comp.Title.FontWeight = 'bold';
    catch
        warning('Impossible de modifier la taille du titre comparative');
    end
    
    % Rotation des étiquettes X
    try
        % Pour les versions récentes de MATLAB
        h_comp.XDisplayLabels = cellfun(@(x) strrep(x, ' ', {' '; ''}), allVars, 'UniformOutput', false);
    catch
        % Méthode alternative si la propriété n'existe pas
        warning('Impossible de modifier les étiquettes X de la heatmap comparative');
    end
    
    exportgraphics(h_comp, fullfile(outdir, "07_Heatmap_Comparative_ALL_"+ts+".png"), 'Resolution',300);
    
    % Sauvegarder les données comparatives
    comp_table = array2table(all_profiles, 'VariableNames', allVars, 'RowNames', all_labels);
    writetable(addvars(comp_table, (1:size(all_profiles,1))', 'Before', 1, 'NewVariableNames', 'ProfileID'), ...
               fullfile(outdir, "08_Comparative_Profiles_"+ts+".csv"), ...
               'WriteRowNames', true);
end

disp("✅ Clustering global + par condition + heatmaps sauvegardés dans : " + outdir);

%% ===================== FONCTIONS LOCALES =====================

function sumd = elbow_sumd(pcs, kmax)
% Retourne la somme des distances intra-cluster pour k=1..kmax
kmax = max(1, kmax);
sumd = zeros(1,kmax);
rng(42);
for ktest = 1:kmax
    if size(pcs,1) >= ktest
        [~,~, dist] = kmeans(pcs, ktest, 'Replicates',10, 'Display','off');
        sumd(ktest) = sum(dist);
    else
        sumd(ktest) = NaN; % pas assez d'observations
    end
end
% Remplace NaN par interpolation simple si besoin
if any(isnan(sumd))
    valid = ~isnan(sumd);
    sumd(~valid) = interp1(find(valid), sumd(valid), find(~valid), 'linear', 'extrap');
end
end

function [k_auto, dist_line] = elbow_k(sumd)
% Trouve le "coude" via distance max à la corde
x = 1:numel(sumd); y = sumd(:)';
p1 = [x(1) y(1)]; p2 = [x(end) y(end)];
num = abs((p2(2)-p1(2))*x - (p2(1)-p1(1))*y + p2(1)*p1(2) - p2(2)*p1(1));
den = sqrt( (p2(2)-p1(2))^2 + (p2(1)-p1(1))^2 );
dist_line = num ./ max(den, eps);
[~, k_auto_rel] = max(dist_line(2:end-1));
k_auto = k_auto_rel + 1;
k_auto = max(2, k_auto); % évite k=1 par défaut
end

function AgeStats = annotate_age_on_clusters(ax, pcs, idxCluster, Cc, agesMonths)
% Affiche pour chaque cluster : "μ ± σ ans" près du centroïde, avec halo.
% Retourne aussi un tableau AgeStats (Cluster, N, AgeMonths_mean, AgeMonths_sd, AgeYears_label)
k = size(Cc,1);
AgeStats = table('Size',[k 5], ...
    'VariableTypes', {'double','double','double','double','string'}, ...
    'VariableNames', {'Cluster','N','AgeMonths_mean','AgeMonths_sd','AgeYears_label'});
AgeStats.Cluster = (1:k)';

hold(ax, 'on');
dx = 0.12; dy = 0.12;   % petit décalage
halo = 0.05;            % rayon du halo

for i = 1:k
    idx = (idxCluster == i);
    ages_m = agesMonths(idx);
    n_i  = numel(ages_m);
    mu_m = mean(ages_m, 'omitnan');
    sd_m = std(ages_m, 'omitnan');

    AgeStats.N(i)              = n_i;
    AgeStats.AgeMonths_mean(i) = mu_m;
    AgeStats.AgeMonths_sd(i)   = sd_m;

    mu_y = mu_m/12; sd_y = sd_m/12;
    lbl = sprintf('%.1f \xB1 %.1f ans', mu_y, sd_y); % \xB1 = ±
    AgeStats.AgeYears_label(i) = lbl;

    if all(isfinite(Cc(i,:)))
        x = Cc(i,1) + dx;
        y = Cc(i,2) + dy;
        txt = "\leftarrow " + string(lbl);

        % Halo blanc (8 copies autour)
        for ang = linspace(0, 2*pi, 8)
            text(ax, x + halo*cos(ang), y + halo*sin(ang), txt, ...
                'HorizontalAlignment','left','VerticalAlignment','bottom', ...
                'FontSize',12,'FontWeight','bold','Color','w','Clipping','on');
        end
        % Texte principal
        text(ax, x, y, txt, ...
            'HorizontalAlignment','left','VerticalAlignment','bottom', ...
            'FontSize',12,'FontWeight','bold','Color','k','Clipping','on');
    end
end
end

function grp = age_group_from_months(m)
% Retourne la catégorie d'âge en fonction de l'âge en mois.
% Bornes utilisées : <72 = JeunesEnfants ; 72-144 = Enfants ;
% 144-216 = Adolescents ; >216 = Adultes
if isnan(m)
    grp = "NA";
elseif m <= 72
    grp = "JeunesEnfants";
elseif m <= 144
    grp = "Enfants";
elseif m <= 216
    grp = "Adolescents";
else
    grp = "Adultes";
end
end

function print_cluster_composition(meta_tbl, idxCluster, cond_list)
% Affiche dans la console, pour chaque cluster :
% - Nombre d'observations par groupe d'âge
% - (Si plusieurs conditions dans meta_tbl) nb de conditions distinctes + répartition par condition

clusters = unique(idxCluster(:))';
ageCats = ["JeunesEnfants","Enfants","Adolescents","Adultes"];

hasMultipleConds = numel(unique(meta_tbl.Condition)) > 1;

for c = clusters
    idx = (idxCluster == c);
    sub = meta_tbl(idx, :);

    % Comptage par groupe d'âge
    countsAge = zeros(size(ageCats));
    for a = 1:numel(ageCats)
        countsAge(a) = sum(sub.AgeGroup == ageCats(a));
    end

    fprintf('Cluster %d:\n', c);
    fprintf('  Par groupe d''âge : ');
    for a = 1:numel(ageCats)
        fprintf('%s=%d%s', ageCats(a), countsAge(a), iff(a<numel(ageCats), ', ', ''));
    end
    fprintf('\n');

    % Conditions
    if hasMultipleConds
        [uConds, ~, ic] = unique(sub.Condition);
        nDistinct = numel(uConds);
        fprintf('  Conditions représentées : %d (', nDistinct);
        % breakdown
        for k = 1:numel(uConds)
            cnt = sum(ic==k);
            fprintf('%s=%d%s', uConds(k), cnt, iff(k<numel(uConds), ', ', ''));
        end
        fprintf(')\n');
    end
end
end

function s = iff(cond, a, b)
% petit utilitaire inline if
if cond, s = a; else, s = b; end
end

function cmap = blueyellow(n)
% Crée une colormap bleu-blanc-jaune centrée sur 0
if nargin < 1, n = 256; end

% Créer le dégradé : Bleu foncé -> Blanc -> Jaune
% Première moitié : Bleu vers blanc
r1 = linspace(0, 1, n/2)';      % Rouge: de 0 à 1
g1 = linspace(0, 1, n/2)';      % Vert: de 0 à 1  
b1 = ones(n/2, 1);              % Bleu: reste à 1

% Deuxième moitié : Blanc vers jaune
r2 = ones(n/2, 1);              % Rouge: reste à 1
g2 = ones(n/2, 1);              % Vert: reste à 1
b2 = linspace(1, 0, n/2)';      % Bleu: de 1 à 0

% Assembler la colormap complète
cmap = [r1, g1, b1; r2, g2, b2];
end