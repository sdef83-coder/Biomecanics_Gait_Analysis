%% CLUSTERING SPATIO-TEMPOREL — PUBLICATION-GRADE (k-means robuste, nPC >= 70% var.)
% - ACP commune (z-score), clustering sur nPC_final (variance cumulée >= 70%)
% - Viz PC1–PC2 (gradient d'âge), mais partition apprise sur nPC_final
% - Sélection de k par vote multi-critères (Sil, CH, DBI, Gap-1SE) + Stabilité (bootstrap ARI)
% - Rapports "littéraires" : z & unités réelles, Cohen's d, tests (ANOVA/Kruskal) + FDR
% - Analyses GLOBAL + PAR CONDITION (Plat/Medium/High), exports CSV/PNG

clc; clear; close all;

% === Dossier d'E/S ===
root_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result';
cd(root_path)
load('SpatioTemporalDATA.mat') % attend SpatioTemporalDATA.ALL.(Plat|Medium|High)

outdir = fullfile(root_path, 'Fig', 'Clustering');
if ~exist(outdir, 'dir'); mkdir(outdir); end
ts = string(datetime('now','Format','yyyyMMdd_HHmm')); % horodatage

% === Domaines de clustering (5 catégories) ===
domains = struct();

domains.Pace = {
    'vitFoulee',      'Gait speed (m.s^{-1})',     'Mean';
    'NormCadence',    'Norm Cadence (ua)',         'Mean';
    'NormStepLength', 'Norm Step length (ua)',     'Mean';
    'NormWalkRatio',  'Norm WR (ua)',              'Mean'
};

domains.Rhythm = {
    'DoubleSupport',  'Double support time (%)',   'Mean'
};

domains.Stability = {
    'LargeurPas',     'Stride width (cm)',         'Mean';
    'MoS_ML_HS_pL0',  'MoS ML HS (%L0)',           'Mean';
    'MoS_AP_HS_pL0',  'MoS AP HS (%L0)',           'Mean'
};

domains.Asymmetry = {
    'distFoulee',     'Stride length (m)',         'SI';
    'DoubleSupport',  'Double support time (%)',   'SI';
    'LargeurPas',     'Stride width (cm)',         'SI';
};

domains.Variability = {
    'distFoulee',     'Stride length (m)',         'CV';
    'LargeurPas',     'Stride width (cm)',         'CV';
    'vitFoulee',      'Gait speed (m.s^{-1})',     'CV'  
};

domains.Smoothness = {
   'STERN_SPARC_Magnitude','STERN SPARC Magnitude (ua)','Mean'
};

% === Mapping "nom technique" → "nom lisible" (ce que connaît la table) ===
nameMap = containers.Map;
nameMap('vitFoulee')               = 'Gait speed (m.s^{-1})';
nameMap('distFoulee')              = 'Stride length (m)';
nameMap('NormStepLength')          = 'Norm Step length (ua)';
nameMap('tempsFoulee')             = 'Stride time (s)';
nameMap('vitCadencePasParMinute')  = 'Cadence (step.min^{-1})';
nameMap('NormCadence')             = 'Norm Cadence (ua)';
nameMap('NormWalkRatio')           = 'Norm WR (ua)';
nameMap('LargeurPas')              = 'Stride width (cm)';
nameMap('DoubleSupport')           = 'Double support time (%)';
nameMap('MoS_AP_HS_pL0')           = 'MoS AP HS (%L0)';
nameMap('MoS_ML_HS_pL0')           = 'MoS ML HS (%L0)';
nameMap('COM_SPARC_Magnitude')     = 'COM SPARC Magnitude (ua)';
nameMap('COM_LDLJ_Magnitude')      = 'COM LDLJ Magnitude (ua)';
nameMap('STERN_SPARC_Magnitude')   = 'STERN SPARC Magnitude (ua)';
nameMap('STERN_LDLJ_Magnitude')    = 'STERN LDLJ Magnitude (ua)';

% === Construction automatique de allVars depuis les domaines ===
allVars = {};
domainNames = fieldnames(domains);

for d = 1:numel(domainNames)
    dVars = domains.(domainNames{d});
    for v = 1:size(dVars, 1)
        varTech = dVars{v,1};  % ex: 'vitFoulee'
        varType = dVars{v,3};  % 'Mean', 'SI', ou 'CV'
        
        % Construire le nom complet : "Mean_Gait speed (m.s^{-1})"
        baseName = nameMap(varTech);
        fullVarName = [varType '_' baseName];
        
        % Éviter les doublons
        if ~ismember(fullVarName, allVars)
            allVars{end+1} = fullVarName; %#ok<SAGROW>
        end
    end
end

fprintf('\n=== VARIABLES SÉLECTIONNÉES POUR LE CLUSTERING ===\n');
fprintf('Total : %d variables réparties en 5 domaines\n', numel(allVars));
for d = 1:numel(domainNames)
    dName = domainNames{d};
    nVars = size(domains.(dName), 1);
    fprintf('  - %s : %d variables\n', dName, nVars);
end
fprintf('\nListe complète :\n');
disp(allVars');

% === Version "affichage" des noms (sans unités) ===
allVars_disp = allVars;
for i = 1:numel(allVars_disp)
    % cas spécifique : CV_Gait speed (m.s^{-1}) -> CV_Stride speed
    if strcmp(allVars_disp{i}, 'CV_Gait speed (m.s^{-1})')
        allVars_disp{i} = 'CV_Stride speed';
        continue;
    end
    % supprimer les unités entre parenthèses partout ailleurs
    allVars_disp{i} = regexprep(allVars_disp{i}, '\s*\([^)]*\)', '');
    allVars_disp{i} = strtrim(regexprep(allVars_disp{i}, '\s+', ' '));
end

% === Conditions ===
conds = {'Plat','Medium','High'};

% === Concaténation des données ===
dataMat = [];
meta = table();  % Condition, Participant, AgeMonths

for iC = 1:numel(conds)
    cond = conds{iC};
    T = SpatioTemporalDATA.ALL.(cond);
    mask = all(~ismissing(T(:, allVars)), 2);
    Tsel = T(mask, :);
    X = Tsel{:, allVars};
    dataMat = [dataMat; X]; %#ok<AGROW>
    n = size(X,1);
    meta = [meta; table(repmat(string(cond),n,1), Tsel.Participant, Tsel.AgeMonths, ...
        'VariableNames',{'Condition','Participant','AgeMonths'})]; %#ok<AGROW>
end

% Groupes d'âge (catégories fixes)
meta.AgeGroup = arrayfun(@(m) age_group_from_months(m), meta.AgeMonths, 'UniformOutput', false);
meta.AgeGroup = string(meta.AgeGroup);

% === Normalisation & ACP commune ===
dataNorm = zscore(dataMat);
[coeff, score, ~, ~, explained] = pca(dataNorm);   % base commune

% Nombre de PCs à afficher (par exemple les 10 premières ou toutes si moins)
nPC_display = min(10, size(coeff, 2));

% === COURBE DE VARIANCE EXPLIQUÉE ===
fig_var = figure('Color','w');
cumVar = cumsum(explained);

yyaxis left
bar(explained, 'FaceColor', [0.2 0.4 0.8]); hold on;
ylabel('Variance expliquée (%)');

yyaxis right
plot(cumVar, '-o', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.5);
ylabel('Variance cumulée (%)');

xlabel('Composantes principales (PC)');
title('Variance expliquée et cumulée — ACP globale');
grid on;

% Ligne du seuil 70 %
yline(70, '--k', '70%');

exportgraphics(fig_var, fullfile(outdir, "01_VarianceExpliquee_GLOBAL_"+ts+".png"), 'Resolution',300);

% === CERCLE DE CORRÉLATION (PC1-PC2) ===
fig_corr = figure('Color','w'); hold on; axis equal;
th = linspace(0, 2*pi, 100);
plot(cos(th), sin(th), 'k--'); % cercle unité
xlabel(sprintf('PC1 (%.1f%%)', explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', explained(2)));
title('Cercle de corrélation (variables originales vs PC1-PC2)');
grid on;

% Corrélations variables vs composantes
corr_vars = corr(dataNorm, score(:,1:2));

% Tracés des flèches
for i = 1:size(corr_vars,1)
    quiver(0,0, corr_vars(i,1), corr_vars(i,2), 0, 'LineWidth',1.2, 'MaxHeadSize',0.1);
    % ici on peut afficher la version sans unités
    text(corr_vars(i,1)*1.1, corr_vars(i,2)*1.1, allVars_disp{i}, ...
        'FontSize',9,'Interpreter','none');
end
xlim([-1.1 1.1]); ylim([-1.1 1.1]);

exportgraphics(fig_corr, fullfile(outdir, "02_CorrelationCircle_PC1PC2_"+ts+".png"), 'Resolution',300);

% === CHOIX DU NOMBRE DE PC POUR LE CLUSTERING (>= 70% de variance) ===
var_threshold = 70; %
cumVar = cumsum(explained);
nPC_final = find(cumVar >= var_threshold, 1, 'first');
if isempty(nPC_final), nPC_final = min(3, size(score,2)); end
nPC_final = max(2, nPC_final); % minimum 2 pour garder une structure
fprintf('nPC_final choisi: %d PCs (variance cumulée = %.1f%%)\n', nPC_final, cumVar(nPC_final));

% === HEATMAP DES LOADINGS ===
% Matrice des loadings (corrélations entre variables et PCs sélectionnées)
loadings = corr(dataNorm, score(:, 1:nPC_final));

fig_heatmap_loadings = figure('Color','w', 'Position', [100 100 1200 800]);
h_load = heatmap(arrayfun(@(i) sprintf('PC%d (%.1f%%)', i, explained(i)), 1:nPC_final, 'UniformOutput', false), ...
                 allVars_disp, loadings, ...
                 'Title', sprintf('Loadings des variables sur les %d PCs sélectionnées (≥70%% variance)', nPC_final), ...
                 'XLabel', 'Composantes principales', ...
                 'YLabel', 'Variables spatio-temporelles', ...
                 'Colormap', redblue(256), ...
                 'ColorLimits', [-1, 1]);
h_load.FontSize = 9;
h_load.CellLabelFormat = '%.2f';
h_load.GridVisible = 'on';

exportgraphics(fig_heatmap_loadings, fullfile(outdir, "01b_Heatmap_Loadings_SelectedPCs_"+ts+".png"), 'Resolution',300);

% === PAIR PLOTS DES PCs SELECTIONNEES ===
nPC_pairplot = min(5, nPC_final);

fig_pairplot = figure('Color','w', 'Position', [50 50 1400 1400]);

% Créer une grille de subplots
for i = 1:nPC_pairplot
    for j = 1:nPC_pairplot
        subplot(nPC_pairplot, nPC_pairplot, (i-1)*nPC_pairplot + j);
        
        if i == j
            % Diagonale : histogramme de la PC
            histogram(score(:,i), 30, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none');
            title(sprintf('PC%d (%.1f%%)', i, explained(i)), 'FontSize', 10);
            ylabel('Fréquence');
            grid on;
        else
            % Hors diagonale : scatter plot PC_j vs PC_i avec gradient d'âge
            scatter(score(:,j), score(:,i), 25, meta.AgeMonths/12, 'filled', 'MarkerFaceAlpha', 0.6);
            colormap(gca, parula);
            xlabel(sprintf('PC%d (%.1f%%)', j, explained(j)), 'FontSize', 9);
            ylabel(sprintf('PC%d (%.1f%%)', i, explained(i)), 'FontSize', 9);
            grid on;
            
            % Colorbar seulement pour le dernier subplot de chaque ligne
            if j == nPC_pairplot
                cb = colorbar('eastoutside');
                cb.Label.String = 'Âge (ans)';
                cb.FontSize = 8;
            end
        end
    end
end

sgtitle(sprintf('Pairplot des %d premières PCs avec gradient d''âge', nPC_pairplot), ...
        'FontSize', 14, 'FontWeight', 'bold');

exportgraphics(fig_pairplot, fullfile(outdir, "02b_Pairplot_PCs_AgeGradient_"+ts+".png"), 'Resolution',300);

% Pour la visualisation, on conserve PC1–PC2 ; pour le clustering, on utilise 1:nPC_final
Xpcs_all_viz = score(:,1:2);
Xpcs_all     = score(:,1:nPC_final);

%% ========== (A) CLUSTERING PAR K-MEANS ANALYSE GLOBALE (toutes conditions) ==========
rng(42);

% --- Sélection robuste de k (vote multi-critères + bootstrap ARI) ---
kRange = 2:10;
[sil_vals,ch_vals,db_vals,gap_vals,gap_th] = criteria_curves(Xpcs_all,kRange);
ari_mean = bootstrap_ari(Xpcs_all,kRange,100);

sil_n = rescale(sil_vals,0,1);
ch_n  = rescale(ch_vals ,0,1);
db_n  = 1 - rescale(db_vals,0,1);
gap_n = rescale(gap_vals,0,1);
ari_n = rescale(ari_mean,0,1);
score_global = 0.30*sil_n + 0.20*ch_n + 0.20*db_n + 0.10*gap_n + 0.20*ari_n;
[~,best_idx] = max(score_global);
k_global = kRange(best_idx);

% --- Courbe du coude (WCSS) ---
WCSS_vals = zeros(size(kRange));
for kk = 1:numel(kRange)
    k = kRange(kk);
    [idx_tmp, Ctmp] = kmeans(Xpcs_all, k, ...
        'Replicates', 20, 'MaxIter', 300, 'Display', 'off');

    % distance de chaque point à SON centroïde
    d2 = sum((Xpcs_all - Ctmp(idx_tmp, :)).^2, 2);  % n×1
    WCSS_vals(kk) = sum(d2);
end
fig_elbow = figure('Color','w');
plot(kRange, WCSS_vals, '-o', 'LineWidth',1.5);
xlabel('Nombre de clusters (k)');
ylabel('Somme des distances intra-clusters (WCSS)');
title('Méthode du coude — Données globales');
grid on;
exportgraphics(fig_elbow, fullfile(outdir, "03a_ElbowPlot_GLOBAL_"+ts+".png"), 'Resolution',300);

% Figures critères
fig_kcrit = figure('Color','w');
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(kRange,sil_vals,'-o'); grid on; title('Silhouette ↑'); xlabel('k');
nexttile; plot(kRange,ch_vals ,'-o'); grid on; title('Calinski–Harabasz ↑'); xlabel('k');
nexttile; plot(kRange,db_vals ,'-o'); grid on; title('Davies–Bouldin ↓'); xlabel('k');
nexttile; plot(kRange,gap_vals,'-o'); grid on; title('Gap ↑ (règle 1-SE)'); xlabel('k'); hold on; yline(gap_th,'k--','1-SE');
nexttile; plot(kRange,ari_mean,'-o'); grid on; title('Stabilité (ARI) ↑'); xlabel('k');
nexttile; plot(kRange,score_global,'-o'); hold on; plot(k_global, score_global(best_idx),'rp','MarkerFaceColor','r');
grid on; title(sprintf('Score global ↑ (k*=%d)',k_global)); xlabel('k');
exportgraphics(fig_kcrit, fullfile(outdir, "03b_Kselection_multiCriteria_GLOBAL_"+ts+".png"), 'Resolution',300);

% === ÉVALUATION COMPARATIVE (k=2→10) — TABLEAU SYNTHÉTIQUE ===
kRange_valid = 2:10;
validity = table('Size',[numel(kRange_valid) 7], ...
    'VariableNames', {'k','CH','Silhouette','DB','Gap','DNg','DNs'}, ...
    'VariableTypes', {'double','double','double','double','double','double','double'});

rng(42);
for ii = 1:numel(kRange_valid)
    k = kRange_valid(ii);
    if size(Xpcs_all,1) < k, continue; end
    idx = kmeans(Xpcs_all, k, 'Replicates',40,'MaxIter',400,'Display','off');
    try SilMean = mean(silhouette(Xpcs_all, idx)); catch, SilMean = NaN; end
    EC_ch  = evalclusters(Xpcs_all,'kmeans','CalinskiHarabasz','KList',k);
    EC_db  = evalclusters(Xpcs_all,'kmeans','DaviesBouldin','KList',k);
    EC_gap = evalclusters(Xpcs_all,'kmeans','gap','KList',k);

    % Distances inter/intra
    D = pdist2(Xpcs_all, Xpcs_all);
    intra = mean(arrayfun(@(c) mean(D(idx==c,idx==c),'all','omitnan'), 1:k));
    inter = mean(arrayfun(@(c) mean(D(idx==c,idx~=c),'all','omitnan'), 1:k));

    validity(ii,:) = {k, EC_ch.CriterionValues, SilMean, EC_db.CriterionValues, EC_gap.CriterionValues, inter, intra};
end

disp('=== Pertinence des clusters (2→10) ===');
disp(validity);

% Export CSV uniquement
out_csv = fullfile(outdir, "VALIDITY_INDEXES_GLOBAL_"+ts+".csv");
writetable(validity, out_csv);
fprintf('✅ Export CSV : %s\n', out_csv);

% --- Clustering final (global) sur nPC_final ---
[idxCluster_global, Cc_global_npc] = kmeans(Xpcs_all, k_global, ...
    'Replicates',50,'MaxIter',500,'Display','off');
Cc_global_pc12 = Cc_global_npc(:,1:2);  % centroïdes projetés sur PC1–PC2

% === TRAÇABILITÉ DES INDIVIDUS PAR CLUSTER ===
ClusterTrace = table();
ClusterTrace.Participant = meta.Participant;
ClusterTrace.Condition   = meta.Condition;
ClusterTrace.AgeMonths   = meta.AgeMonths;
ClusterTrace.AgeGroup    = meta.AgeGroup;
ClusterTrace.Cluster     = idxCluster_global;
ClusterTrace.PC1         = Xpcs_all_viz(:,1);
ClusterTrace.PC2         = Xpcs_all_viz(:,2);

% Tri par cluster
ClusterTrace = sortrows(ClusterTrace, "Cluster");

% Export CSV
out_trace = fullfile(outdir, "06_ClusterTraceability_GLOBAL_"+ts+".csv");
writetable(ClusterTrace, out_trace);
fprintf('✅ Export traçabilité participants : %s\n', out_trace);

% Définir les markers par condition
markerMap = containers.Map({'Plat', 'Medium', 'High'}, {'o', 's', '^'});
markerSizes = 60; % Taille des markers

fig_agegrad_shapes = figure('Color','w', 'Position', [100 100 1000 800]); 
hold on;

% Tracer chaque condition avec sa forme spécifique
conds_unique = unique(meta.Condition);
for iC = 1:numel(conds_unique)
    cond = conds_unique(iC);
    idx_cond = meta.Condition == cond;
    
    scatter(Xpcs_all_viz(idx_cond, 1), Xpcs_all_viz(idx_cond, 2), ...
            markerSizes, meta.AgeMonths(idx_cond)/12, ...
            'filled', markerMap(char(cond)), ...
            'MarkerFaceAlpha', 0.7, ...
            'DisplayName', char(cond));
end

colormap(parula); 
cb = colorbar; 
cb.Label.String = 'Âge (années)';
cb.FontSize = 11;

xlabel(sprintf('PC1 (%.1f%%)', explained(1)), 'FontSize', 12); 
ylabel(sprintf('PC2 (%.1f%%)', explained(2)), 'FontSize', 12);
title(sprintf('PC1–PC2 avec gradient d''âge et formes par surface\nClustering sur %d PC(s), k=%d', ...
              nPC_final, k_global), 'FontSize', 13);
grid on;

% Tracer les contours des clusters
cols = lines(k_global);
for i = 1:k_global
    pts2 = Xpcs_all_viz(idxCluster_global==i, :);
    if size(pts2,1) >= 3
        K2 = convhull(pts2(:,1), pts2(:,2));
        plot(pts2(K2,1), pts2(K2,2), '-', 'Color', cols(i,:), ...
             'LineWidth', 2, 'HandleVisibility', 'off');
    end
end

% Centroïdes
plot(Cc_global_pc12(:,1), Cc_global_pc12(:,2), 'kx', ...
     'MarkerSize', 14, 'LineWidth', 3, 'DisplayName', 'Centroïdes');
for i = 1:k_global
    text(Cc_global_pc12(i,1), Cc_global_pc12(i,2), sprintf('  C%d', i), ...
         'FontWeight', 'bold', 'Color', 'k', 'VerticalAlignment', 'middle', ...
         'FontSize', 11);
end

% Légende pour les surfaces
legend('Location', 'best', 'FontSize', 10);

% Annotations d'âge
AgeStats_global = annotate_age_on_clusters(gca, Xpcs_all_viz, idxCluster_global, ...
                                           Cc_global_pc12, meta.AgeMonths);

exportgraphics(fig_agegrad_shapes, ...
               fullfile(outdir, "04_PC1PC2_AgeGradient_Shapes_GLOBAL_k"+k_global+"_"+ts+".png"), ...
               'Resolution', 300);

% --- Profils moyens: z & unités réelles, Cohen's d, tests + FDR ---
[meanZ_global, meanRAW_global, cohenD_global, statsTable_global] = ...
    cluster_profiles_and_stats(dataMat, dataNorm, idxCluster_global, allVars);

% --- Sauvegardes GLOBAL ---
writetable(AgeStats_global, fullfile(outdir, "Age_Stats_GLOBAL_"+ts+".csv"));
writetable(addvars(meanZ_global,(1:k_global)','Before',1,'NewVariableNames','ClusterID'), ...
           fullfile(outdir, "ProfilMoyen_zscore_GLOBAL_"+ts+".csv"), 'WriteRowNames',true);
writetable(addvars(meanRAW_global,(1:k_global)','Before',1,'NewVariableNames','ClusterID'), ...
           fullfile(outdir, "ProfilMoyen_unitesReelles_GLOBAL_"+ts+".csv"), 'WriteRowNames',true);
writetable(addvars(cohenD_global,(1:k_global)','Before',1,'NewVariableNames','ClusterID'), ...
           fullfile(outdir, "CohenD_GLOBAL_"+ts+".csv"), 'WriteRowNames',true);
writetable(statsTable_global, fullfile(outdir,"Stats_omnibus_GLOBAL_"+ts+".csv"));

% --- Heatmap z-score (UNE SEULE FIGURE, annotations intégrées) ---
fig_heat = figure('Color','w');
heatmap_data_global = table2array(meanZ_global);
clim = max(abs(heatmap_data_global(:)));
rowLabels = arrayfun(@(i) ...
    sprintf('C%d (n=%d, %.1f ans)', i, sum(idxCluster_global==i), AgeStats_global.AgeMonths_mean(i)/12), ...
    1:k_global, 'UniformOutput', false);

h = heatmap(allVars_disp, rowLabels, heatmap_data_global, ...
    'Title', sprintf('Profils moyens z-score — Global (k=%d, nPC=%d)', k_global, nPC_final), ...
    'XLabel','Variables spatio-temporelles','YLabel','Clusters (taille, âge)', ...
    'Colormap', blueyellow(256), 'ColorLimits', [-clim, clim]);
h.FontSize = 10; h.CellLabelFormat = '%.2f'; h.GridVisible = 'on';
exportgraphics(fig_heat, fullfile(outdir, "05_Heatmap_Profils_GLOBAL_k"+k_global+"_"+ts+".png"), 'Resolution',300);

% --- Qualité & composition ---
[sil_each, CH_best, DB_best, GAP_best, EntTable_global] = ...
    quality_and_composition(Xpcs_all, idxCluster_global, meta, k_global, ch_vals, db_vals, gap_vals, kRange);
writetable( table((1:numel(sil_each))', sil_each, idxCluster_global, ...
    'VariableNames', {'Obs','Silhouette','Cluster'}) , ...
    fullfile(outdir,"Silhouette_indiv_GLOBAL_"+ts+".csv") );
writetable(EntTable_global, fullfile(outdir,"Qualite_Clusters_GLOBAL_"+ts+".csv"));
fprintf('\n=== QUALITE (GLOBAL) ===\n');
fprintf('nPC=%d | k*=%d | Sil=%.3f | CH=%.1f | DBI=%.3f | GAP=%.3f | ARI_boot=%.3f\n', ...
    nPC_final, k_global, mean(sil_each,'omitnan'), CH_best, DB_best, GAP_best, ari_mean(kRange==k_global));
print_cluster_composition(meta, idxCluster_global, conds);

% === TESTS STATISTIQUES ENTRE CLUSTERS ===
fprintf('\n=== Tests statistiques entre clusters (GLOBAL) ===\n');

cluster_labels = unique(idxCluster_global);
nClust = numel(cluster_labels);
pvals = nan(1, numel(allVars));
cohensD = nan(1, numel(allVars));
testType = strings(numel(allVars),1);

for v = 1:numel(allVars)
    data1 = dataMat(idxCluster_global == cluster_labels(1), v);
    data2 = dataMat(idxCluster_global == cluster_labels(2), v);

    % Test de normalité
    if all(~isnan(data1)) && all(~isnan(data2))
        isNorm = all([kstest((data1-mean(data1))/std(data1)) == 0, ...
                      kstest((data2-mean(data2))/std(data2)) == 0]);
    else
        isNorm = false;
    end

    % Test statistique
    if isNorm
        [~, p] = ttest2(data1, data2);
        testType(v) = "t-test";
    else
        p = ranksum(data1, data2);
        testType(v) = "Mann-Whitney";
    end
    pvals(v) = p;

    % Taille d'effet (Cohen's d)
    m1 = mean(data1, 'omitnan'); m2 = mean(data2, 'omitnan');
    s1 = std(data1, 'omitnan'); s2 = std(data2, 'omitnan');
    s_pooled = sqrt(((numel(data1)-1)*s1^2 + (numel(data2)-1)*s2^2) / ...
                    (numel(data1)+numel(data2)-2));
    cohensD(v) = (m1 - m2) / s_pooled;
end

% Correction Bonferroni
nTests = numel(pvals);
pvals_bonf = pvals * nTests;
pvals_bonf(pvals_bonf > 1) = 1; % borne max à 1

% Indicateurs de significativité
Sig = repmat("ns", numel(pvals), 1);
Sig(pvals_bonf < 0.05) = "*";
Sig(pvals_bonf < 0.01) = "**";
Sig(pvals_bonf < 0.001) = "***";

% Résumé tableau avec Bonferroni + type de test
StatsClusters = table(allVars_disp', testType, pvals', pvals_bonf', cohensD', Sig, ...
    'VariableNames', {'Variable','Test','p_value','p_Bonferroni','Cohens_d','Sig'});
disp(StatsClusters)

% Export CSV
writetable(StatsClusters, fullfile(outdir, "ClusterComparison_STATS_BONFERRONI_"+ts+".csv"));

%% ========== (B) CLUSTERING PAR K-MEANS ANALYSES PAR CONDITION (Plat / Medium / High) ==========
independentPCA = false;  % false = ACP commune (comparabilité), true = ACP par condition

for iC = 1:numel(conds)
    cond = conds{iC};
    idxCond = meta.Condition==cond;
    nCond = sum(idxCond);
    if nCond < 6
        warning('Condition %s: trop peu d''observations (%d). Clustering sauté.', cond, nCond);
        continue;
    end

    if independentPCA
        Xnorm_c = dataNorm(idxCond, :);
        [~, score_c, ~, ~, explained_c] = pca(Xnorm_c);
        pcs_c_viz = score_c(:,1:2);              % viz
        pcs_c     = score_c(:,1:min(nPC_final,size(score_c,2))); % clustering
        expl1 = explained_c(1); expl2 = explained_c(2);
    else
        pcs_c_viz = Xpcs_all_viz(idxCond,:);     % projection 2D commune
        pcs_c     = Xpcs_all(idxCond,:);         % nPC_final commun
        expl1 = explained(1); expl2 = explained(2);
        Xnorm_c = dataNorm(idxCond, :);
    end
    meta_c = meta(idxCond,:);

    % === ÉVALUATION COMPARATIVE (k=2→10) — TABLEAU SYNTHÉTIQUE ===
    kRange_valid_c = 2:min(10, max(2, nCond-1));
    validity_c = table('Size',[numel(kRange_valid_c) 7], ...
        'VariableNames', {'k','CH','Silhouette','DB','Gap','DNg','DNs'}, ...
        'VariableTypes', {'double','double','double','double','double','double','double'});

    rng(42);
    for jj = 1:numel(kRange_valid_c)
        kk = kRange_valid_c(jj);
        if size(pcs_c,1) < kk, continue; end
        idx_tmp = kmeans(pcs_c, kk, 'Replicates',30,'MaxIter',400,'Display','off');
        try SilMean_c = mean(silhouette(pcs_c, idx_tmp)); catch, SilMean_c = NaN; end
        EC_ch_c  = evalclusters(pcs_c,'kmeans','CalinskiHarabasz','KList',kk);
        EC_db_c  = evalclusters(pcs_c,'kmeans','DaviesBouldin','KList',kk);
        EC_gap_c = evalclusters(pcs_c,'kmeans','gap','KList',kk);
        % Distances inter/intra
        D_c = pdist2(pcs_c, pcs_c);
        intra_c = mean(arrayfun(@(cl) mean(D_c(idx_tmp==cl,idx_tmp==cl),'all','omitnan'), 1:kk));
        inter_c = mean(arrayfun(@(cl) mean(D_c(idx_tmp==cl,idx_tmp~=cl),'all','omitnan'), 1:kk));
        validity_c(jj,:) = {kk, EC_ch_c.CriterionValues, SilMean_c, EC_db_c.CriterionValues, EC_gap_c.CriterionValues, inter_c, intra_c};
    end
    writetable(validity_c, fullfile(outdir, "B00_VALIDITY_INDEXES_"+cond+"_"+ts+".csv"));

    % Sélection robuste de k (condition)
    kRange_c = 2:min(10, max(2, nCond-1));
    [sil_v_c,ch_v_c,db_v_c,gap_v_c,gap_th_c] = criteria_curves(pcs_c,kRange_c);
    ari_mean_c = bootstrap_ari(pcs_c,kRange_c,80);

    sil_n = rescale(sil_v_c,0,1);
    ch_n  = rescale(ch_v_c ,0,1);
    db_n  = 1 - rescale(db_v_c,0,1);
    gap_n = rescale(gap_v_c,0,1);
    ari_n = rescale(ari_mean_c,0,1);

    score_c = 0.30*sil_n + 0.20*ch_n + 0.20*db_n + 0.10*gap_n + 0.20*ari_n;
    [~,best_idx_c] = max(score_c);
    k_c = kRange_c(best_idx_c);

    % --- Courbe du coude (WCSS) --- (corrigée pour la condition)
    WCSS_vals_c = zeros(size(kRange_c));
    for kk = 1:numel(kRange_c)
        ktest = kRange_c(kk);
        [idx_tmp, Ctmp] = kmeans(pcs_c, ktest, ...
            'Replicates', 20, 'MaxIter', 300, 'Display', 'off');
        d2 = sum((pcs_c - Ctmp(idx_tmp, :)).^2, 2);
        WCSS_vals_c(kk) = sum(d2);
    end
    fig_elbow_c = figure('Color','w');
    plot(kRange_c, WCSS_vals_c, '-o', 'LineWidth',1.5);
    xlabel('Nombre de clusters (k)');
    ylabel('Somme des distances intra-clusters (WCSS)');
    title(sprintf('Méthode du coude — %s', cond));
    grid on;
    exportgraphics(fig_elbow_c, fullfile(outdir, "B01_ElbowPlot_"+cond+"_"+ts+".png"), 'Resolution',300);

    % Figures critères (condition)
    fig_kcrit_c = figure('Color','w');
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    nexttile; plot(kRange_c,sil_v_c,'-o'); grid on; title(sprintf('%s: Silhouette ↑',cond)); xlabel('k');
    nexttile; plot(kRange_c,ch_v_c ,'-o'); grid on; title(sprintf('%s: CH ↑',cond)); xlabel('k');
    nexttile; plot(kRange_c,db_v_c ,'-o'); grid on; title(sprintf('%s: DBI ↓',cond)); xlabel('k');
    nexttile; plot(kRange_c,gap_v_c,'-o'); grid on; title(sprintf('%s: Gap ↑ (1-SE)',cond)); xlabel('k'); hold on; yline(gap_th_c,'k--','1-SE');
    nexttile; plot(kRange_c,ari_mean_c,'-o'); grid on; title(sprintf('%s: ARI ↑',cond)); xlabel('k');
    nexttile; plot(kRange_c,score_c,'-o'); hold on; plot(k_c, score_c(best_idx_c),'rp','MarkerFaceColor','r');
    grid on; title(sprintf('%s: Score global (k*=%d, nPC=%d)',cond,k_c,size(pcs_c,2))); xlabel('k');
    exportgraphics(fig_kcrit_c, fullfile(outdir, "B02_Kselection_"+cond+"_"+ts+".png"), 'Resolution',300);

    % Clustering final (condition)
    [idxC, CcC_npc] = kmeans(pcs_c, k_c, 'Replicates',40,'MaxIter',400,'Display','off');
    CcC_pc12 = CcC_npc(:,1:2);

    % Scatter PC1–PC2 + gradient âge + contours
    figc = figure('Color','w'); hold on;
    scatter(pcs_c_viz(:,1), pcs_c_viz(:,2), 42, meta_c.AgeMonths/12, 'filled');
    colormap(parula); cb = colorbar; cb.Label.String = 'Âge (années)';
    xlabel(sprintf('PC1 (%.1f%%)', expl1)); ylabel(sprintf('PC2 (%.1f%%)', expl2));
    title(sprintf('%s — PC1–PC2 (viz). Clustering sur %d PC(s), k=%d', cond, size(pcs_c,2), k_c)); grid on;

    colsC = lines(k_c);
    for ic = 1:k_c
        pts = pcs_c_viz(idxC==ic,:);
        if size(pts,1) >= 3
            K = convhull(pts(:,1), pts(:,2));
            plot(pts(K,1), pts(K,2), '-', 'Color', colsC(ic,:), 'LineWidth',1.5);
        end
    end
    plot(CcC_pc12(:,1), CcC_pc12(:,2), 'kx', 'MarkerSize',12, 'LineWidth',2);
    for ic=1:k_c
        text(CcC_pc12(ic,1), CcC_pc12(ic,2), sprintf('  C%d',ic), ...
            'FontWeight','bold','Color','k','VerticalAlignment','middle');
    end

    AgeStats_c = annotate_age_on_clusters(gca, pcs_c_viz, idxC, CcC_pc12, meta_c.AgeMonths);
    exportgraphics(figc, fullfile(outdir, "B03_PC1PC2_AgeGradient_"+cond+"_k"+k_c+"_"+ts+".png"), 'Resolution',300);

    % Profils moyens, Cohen's d, tests omnibus
    [meanZ_c, meanRAW_c, cohenD_c, statsTable_c] = ...
        cluster_profiles_and_stats(dataMat(idxCond,:), Xnorm_c, idxC, allVars);

    % === version affichage sans unités pour la condition ===
    allVars_disp_c = allVars;
    for i = 1:numel(allVars_disp_c)
        if strcmp(allVars_disp_c{i}, 'CV_Gait speed (m.s^{-1})')
            allVars_disp_c{i} = 'CV_Stride speed';
            continue;
        end
        allVars_disp_c{i} = regexprep(allVars_disp_c{i}, '\s*\([^)]*\)', '');
        allVars_disp_c{i} = strtrim(regexprep(allVars_disp_c{i}, '\s+', ' '));
    end

    % Heatmap z-score
    fig_heat_c = figure('Color','w');
    heatmap_data_c = table2array(meanZ_c);
    clim_c = max(abs(heatmap_data_c(:)));
    rowLabels_c = arrayfun(@(ic) ...
        sprintf('C%d (n=%d, %.1f ans)', ic, sum(idxC==ic), AgeStats_c.AgeMonths_mean(ic)/12), ...
        1:k_c, 'UniformOutput', false);
    h_c = heatmap(allVars_disp_c, rowLabels_c, heatmap_data_c, ...
        'Title', sprintf('Profils moyens z-score — %s (k=%d, nPC=%d)', cond, k_c, size(pcs_c,2)), ...
        'XLabel','Variables spatio-temporelles','YLabel','Clusters (taille, âge)', ...
        'Colormap', blueyellow(256), 'ColorLimits', [-clim_c, clim_c]);
    h_c.FontSize = 10; h_c.CellLabelFormat = '%.2f'; h_c.GridVisible = 'on';
    exportgraphics(fig_heat_c, fullfile(outdir, "B04_Heatmap_Profils_"+cond+"_k"+k_c+"_"+ts+".png"), 'Resolution',300);

    % Qualité & composition
    [sil_each_c, CH_best_c, DB_best_c, GAP_best_c, EntTable_c] = ...
        quality_and_composition(pcs_c, idxC, meta_c, k_c, ch_v_c, db_v_c, gap_v_c, kRange_c);
    writetable( table((1:numel(sil_each_c))', sil_each_c, idxC, ...
        'VariableNames', {'Obs','Silhouette','Cluster'}) , ...
        fullfile(outdir,"B05_Silhouette_indiv_"+cond+"_"+ts+".csv") );
    writetable(EntTable_c, fullfile(outdir,"B4_Qualite_Clusters_"+cond+"_"+ts+".csv"));
    
    % Exports CSV
    writetable(AgeStats_c, fullfile(outdir, "B5_Age_Stats_"+cond+"_"+ts+".csv"));
    writetable(addvars(meanZ_c,(1:k_c)','Before',1,'NewVariableNames','ClusterID'), ...
        fullfile(outdir, "B6_ProfilMoyen_zscore_"+cond+"_"+ts+".csv"), 'WriteRowNames',true);
    writetable(addvars(meanRAW_c,(1:k_c)','Before',1,'NewVariableNames','ClusterID'), ...
        fullfile(outdir, "B7_ProfilMoyen_unitesReelles_"+cond+"_"+ts+".csv"), 'WriteRowNames',true);
    writetable(addvars(cohenD_c,(1:k_c)','Before',1,'NewVariableNames','ClusterID'), ...
        fullfile(outdir, "B8_CohenD_"+cond+"_"+ts+".csv"), 'WriteRowNames',true);
    writetable(statsTable_c, fullfile(outdir,"B9_Stats_omnibus_"+cond+"_"+ts+".csv"));

    % Tests statistiques entre clusters (seulement si k=2)
    if k_c == 2
        fprintf('\n=== Tests statistiques entre clusters — %s (k=2) ===\n', cond);
        cLabs = unique(idxC);
        pvals_c = nan(1, numel(allVars));
        d_cohen_c = nan(1, numel(allVars));
        testType_c = strings(numel(allVars),1);
        for v = 1:numel(allVars)
            x1 = dataMat(idxCond, v); x1 = x1(idxC == cLabs(1));
            x2 = dataMat(idxCond, v); x2 = x2(idxC == cLabs(2));
            ok1 = ~isnan(x1); ok2 = ~isnan(x2); x1 = x1(ok1); x2 = x2(ok2);
            if numel(x1) >= 3 && numel(x2) >= 3
                norm1 = (lillietest((x1-mean(x1))/std(x1)) == 0);
                norm2 = (lillietest((x2-mean(x2))/std(x2)) == 0);
                if norm1 && norm2
                    [~, p] = ttest2(x1, x2);
                    testType_c(v) = "t-test";
                else
                    p = ranksum(x1, x2);
                    testType_c(v) = "Mann-Whitney";
                end
                pvals_c(v) = p;
                m1 = mean(x1,'omitnan'); m2 = mean(x2,'omitnan');
                s1 = std(x1,'omitnan'); s2 = std(x2,'omitnan');
                sp = sqrt(((numel(x1)-1)*s1^2 + (numel(x2)-1)*s2^2) / (numel(x1)+numel(x2)-2));
                d_cohen_c(v) = (m1 - m2) / max(sp, eps);
            end
        end
        nT = sum(~isnan(pvals_c));
        p_bonf_c = pvals_c * max(nT,1);
        p_bonf_c(p_bonf_c > 1) = 1;
        Sig_c = repmat("ns", numel(pvals_c), 1);
        Sig_c(p_bonf_c < 0.05) = "*";
        Sig_c(p_bonf_c < 0.01) = "**";
        Sig_c(p_bonf_c < 0.001) = "***";
        StatsClusters_c = table(allVars_disp_c', testType_c, pvals_c', p_bonf_c', d_cohen_c', Sig_c, ...
            'VariableNames', {'Variable','Test','p_value','p_Bonferroni','Cohens_d','Sig'});
        disp(StatsClusters_c);
        writetable(StatsClusters_c, fullfile(outdir, "B9b_ClusterComparison_"+cond+"_BONFERRONI_"+ts+".csv"));
    end

    % Résumé console
    fprintf('\n=== QUALITE — %s ===\n', cond);
    fprintf('nPC=%d | k*=%d | Sil=%.3f | CH=%.1f | DBI=%.3f | GAP=%.3f | ARI_boot=%.3f\n', ...
        size(pcs_c,2), k_c, mean(sil_each_c,'omitnan'), CH_best_c, DB_best_c, GAP_best_c, ari_mean_c(kRange_c==k_c));
    print_cluster_composition(meta_c, idxC, {cond});

    Trace_c = table(meta_c.Participant, meta_c.Condition, meta_c.AgeMonths, meta_c.AgeGroup, idxC, pcs_c_viz(:,1), pcs_c_viz(:,2), ...
        'VariableNames', {'Participant','Condition','AgeMonths','AgeGroup','Cluster','PC1','PC2'});
    writetable(Trace_c, fullfile(outdir, "B10_ClusterTraceability_"+cond+"_"+ts+".csv"));
end

disp("✅ Analyses GLOBAL + PAR CONDITION terminées. Exports : " + outdir);

%% ===================== FONCTIONS LOCALES =====================

function grp = age_group_from_months(m)
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

function [sil_vals,ch_vals,db_vals,gap_vals,gap_th] = criteria_curves(X,kRange)
sil_vals = zeros(numel(kRange),1);
ch_vals  = zeros(numel(kRange),1);
db_vals  = zeros(numel(kRange),1);
gap_vals = zeros(numel(kRange),1);
for ii = 1:numel(kRange)
    k = kRange(ii);
    idx_tmp = kmeans(X,k,'Replicates',10,'MaxIter',300,'Display','off');
    try
        sil_vals(ii) = mean(silhouette(X, idx_tmp, 'sqeuclidean'));
    catch
        sil_vals(ii) = NaN;
    end
    EC_ch  = evalclusters(X,'kmeans','CalinskiHarabasz','KList',k);
    EC_db  = evalclusters(X,'kmeans','DaviesBouldin'   ,'KList',k);
    EC_gap = evalclusters(X,'kmeans','gap'              ,'KList',k);
    ch_vals(ii)  = EC_ch.CriterionValues;
    db_vals(ii)  = EC_db.CriterionValues;
    gap_vals(ii) = EC_gap.CriterionValues;
end
gap_se = std(gap_vals,'omitnan');
gap_th = max(gap_vals,[],'omitnan') - gap_se;
end

function ari_mean = bootstrap_ari(X,kRange,B)
ari_mean = zeros(numel(kRange),1);
for ii=1:numel(kRange)
    k = kRange(ii);
    idx_ref = kmeans(X,k,'Replicates',50,'MaxIter',500,'Display','off');
    ARI = zeros(B,1);
    for b=1:B
        boot_idx = randsample(size(X,1), size(X,1), true);
        Xb = X(boot_idx,:);
        idx_b = kmeans(Xb,k,'Replicates',20,'MaxIter',300,'Display','off');
        ARI(b) = adjustedRandIndex(idx_ref(boot_idx), idx_b);
    end
    ari_mean(ii) = mean(ARI,'omitnan');
end
end

function [meanZ, meanRAW, cohenD, statsTable] = cluster_profiles_and_stats(dataRAW, dataZ, idx, allVars)
k = max(idx);
meanZ = array2table(zeros(k,numel(allVars)),'VariableNames',allVars,...
    'RowNames', "Cluster_"+string(1:k));
for i=1:k
    meanZ{i,:} = mean(dataZ(idx==i,:),1,'omitnan');
end
mu = mean(dataRAW,1,'omitnan'); sd = std(dataRAW,0,1,'omitnan');
meanRAW = meanZ;
for v=1:numel(allVars)
    meanRAW{:,v} = meanZ{:,v}.*sd(v) + mu(v);
end
cohenD = array2table(zeros(k,numel(allVars)),'VariableNames',allVars,...
    'RowNames', "Cluster_"+string(1:k));
for i=1;k
    Xc = dataRAW(idx==i,:);
    for v=1:numel(allVars)
        cohenD{i,v} = (mean(Xc(:,v),'omitnan') - mu(v)) / (sd(v)+eps);
    end
end
pvals = nan(1,numel(allVars));
for v=1:numel(allVars)
    y = dataRAW(:,v); g = idx;
    ok = ~isnan(y) & ~isnan(g);
    y = y(ok); g = g(ok);
    if numel(y) > 5 && numel(unique(g)) > 1
        ystd = (y - mean(y,'omitnan'))./std(y,[],'omitnan');
        try normok = (lillietest(ystd) == 0); catch, normok = false; end
        if normok
            pvals(v) = anova1(y,g,'off');
        else
            pvals(v) = kruskalwallis(y,g,'off');
        end
    end
end
[~,~,~,qvals] = fdr_bh(pvals, 0.05);
statsTable = table(allVars(:), pvals(:), qvals(:), ...
    'VariableNames', {'Variable','p_omnibus','q_FDR'});
end

function [sil_each, CH_best, DB_best, GAP_best, EntTable] = quality_and_composition(Xpcs, idx, meta_tbl, k, ch_vals, db_vals, gap_vals, kRange)
try
    sil_each = silhouette(Xpcs, idx);
catch
    sil_each = nan(size(idx));
end
CH_best  = ch_vals( kRange==k );
DB_best  = db_vals( kRange==k );
GAP_best = gap_vals(kRange==k );

ageCats = ["JeunesEnfants","Enfants","Adolescents","Adultes"];
ent_age = zeros(k,1);
ent_surf= zeros(k,1);
for c=1:k
    sub = meta_tbl(idx==c,:);
    countsA = arrayfun(@(a) sum(sub.AgeGroup==a), ageCats);
    pA = countsA/sum(countsA); pA = pA(pA>0);
    ent_age(c) = -sum(pA.*log2(pA));
    catsS = categorical(sub.Condition);
    countsS = countcats(catsS);
    pS = countsS/sum(countsS); pS = pS(pS>0);
    ent_surf(c) = -sum(pS.*log2(pS));
end
sil_byC = splitapply(@(x) mean(x,'omitnan'), sil_each, findgroups(idx));
EntTable = table((1:k)', sil_byC, ent_age, ent_surf, ...
    'VariableNames', {'Cluster','SilhouetteMean','EntropyAge','EntropySurface'});
end

function ARI = adjustedRandIndex(labelsTrue, labelsPred)
labelsTrue = labelsTrue(:); labelsPred = labelsPred(:);
n = numel(labelsTrue);
[~,~,l1] = unique(labelsTrue);
[~,~,l2] = unique(labelsPred);
k1 = max(l1); k2 = max(l2);
N = accumarray([l1 l2],1,[k1 k2]);
ni = sum(N,2); nj = sum(N,1);
sumComb = @(x) sum(x.*(x-1)/2);
t1 = sumComb(N(:));
t2 = sumComb(ni);
t3 = sumComb(nj);
t4 = n*(n-1)/2;
exp_index = (t2*t3)/t4;
max_index = (t2 + t3)/2;
ARI = (t1 - exp_index) / (max_index - exp_index + eps);
end

function [h, crit_p, adj_p, q] = fdr_bh(pvals,qlevel)
if nargin<2, qlevel = 0.05; end
p = pvals(:);
[ps, idx] = sort(p);
m = numel(p);
thresh = (1:m)'/m*qlevel;
is_sig = ps <= thresh;
k = find(is_sig,1,'last');
h = false(size(p));
if ~isempty(k), h(idx(1:k)) = true; end
if ~isempty(k), crit_p = ps(k); else, crit_p = NaN; end
adj_p = nan(size(p));
adj_p(idx) = ps .* m ./ (1:m)';
q = adj_p;
end

function AgeStats = annotate_age_on_clusters(ax, pcs, idxCluster, Cc, agesMonths)
k = size(Cc,1);
AgeStats = table('Size',[k 5], ...
    'VariableTypes', {'double','double','double','double','string'}, ...
    'VariableNames', {'Cluster','N','AgeMonths_mean','AgeMonths_sd','AgeYears_label'});
AgeStats.Cluster = (1:k)';

hold(ax, 'on');
dx = 0.12; dy = 0.12; halo = 0.05;

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
    lbl = sprintf('%.1f %s %.1f ans', mu_y, char(177), sd_y);

    if all(isfinite(Cc(i,:)))
        x = Cc(i,1) + dx; y = Cc(i,2) + dy; txt = "← " + string(lbl);
        for ang = linspace(0, 2*pi, 8)
            text(ax, x + halo*cos(ang), y + halo*sin(ang), txt, ...
                'HorizontalAlignment','left','VerticalAlignment','bottom', ...
                'FontSize',12,'FontWeight','bold','Color','w','Clipping','on');
        end
        text(ax, x, y, txt, 'HorizontalAlignment','left','VerticalAlignment','bottom', ...
            'FontSize',12,'FontWeight','bold','Color','k','Clipping','on');
    end
end
end

function print_cluster_composition(meta_tbl, idxCluster, varargin)
clusters = unique(idxCluster(:))';
ageCats = ["JeunesEnfants","Enfants","Adolescents","Adultes"];
hasMultipleConds = numel(unique(meta_tbl.Condition)) > 1;

for c = clusters
    idx = (idxCluster == c);
    sub = meta_tbl(idx, :);
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
    if hasMultipleConds
        [uConds, ~, ic] = unique(sub.Condition);
        fprintf('  Conditions représentées : %d (', numel(uConds));
        for k = 1:numel(uConds)
            cnt = sum(ic==k);
            fprintf('%s=%d%s', uConds(k), cnt, iff(k<numel(uConds), ', ', ''));
        end
        fprintf(')\n');
    end
end
end

function s = iff(cond, a, b)
if cond, s = a; else, s = b; end
end

function cmap = blueyellow(n)
if nargin < 1, n = 256; end
r1 = linspace(0, 1, n/2)'; g1 = linspace(0, 1, n/2)'; b1 = ones(n/2, 1);
r2 = ones(n/2, 1); g2 = ones(n/2, 1); b2 = linspace(1, 0, n/2)';
cmap = [r1, g1, b1; r2, g2, b2];
end

function cmap = redblue(n)
    % Colormap rouge-blanc-bleu pour loadings
    if nargin < 1, n = 256; end
    
    % Rouge -> Blanc
    r1 = ones(n/2, 1);
    g1 = linspace(0, 1, n/2)';
    b1 = linspace(0, 1, n/2)';
    
    % Blanc -> Bleu
    r2 = linspace(1, 0, n/2)';
    g2 = linspace(1, 0, n/2)';
    b2 = ones(n/2, 1);
    
    cmap = [r1, g1, b1; r2, g2, b2];
end