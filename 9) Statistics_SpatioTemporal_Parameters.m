%% ANALYSES STATISTIQUES AUTOMATISÉES - MODÈLES LINÉAIRES MIXTES (LMM)
% Design: Surface (mesures répétées) × Groupe d'âge (inter-sujets)
% Objectifs de correction:
% 1) Interaction: obtenir F/DF/p (pas de NaN) -> test global via coefTest
% 2) Reporter REML vs ML + DFMethod -> paramètres fixés et sauvegardés
% 3) Garder random intercept seul -> (1|Participant)
% 4) Clarifier FDR -> FDR appliquée séparément par "famille" d'effets (Surface / AgeGroup / Interaction)

clc; clear; close all;

%% =============== Chargement des données ===============
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result');
load('SpatioTemporalDATA.mat');   % doit contenir SpatioTemporalDATA

% Dossier de sortie
stats_path = fullfile(pwd, 'Statistical_Analysis_LMM');
if ~exist(stats_path, 'dir'), mkdir(stats_path); end

%% =============== Configuration ===============
surfaces = {'Plat', 'Medium', 'High'};
groups   = {'JeunesEnfants', 'Enfants', 'Adolescents', 'Adultes'};

% Choix méthodologiques
FIT_METHOD = 'REML';            % 'REML' ou 'ML'
DF_METHOD  = 'Satterthwaite';   % si non supporté back off vers 'Residual'

variables_to_test = {
    % --- Paramètres spatio-temporelles ---
    'Mean_Single support time (%)'
    'Mean_Double support time (%)'
    'Mean_BaseOfSupport (cm)'
    'Mean_StepWidth (cm)'
    'Mean_Gait speed (m.s^{-1})'
    'Mean_Stride length (m)'
    'Mean_Stride time (s)'
    'Mean_WalkRatio'
    'Mean_Norm WR (ua)'
    'Mean_Cadence (step.min^{-1})'
    'Mean_Norm Step length (ua)'
    'Mean_Norm Cadence (ua)'
    'Mean_Norm StepWidth (ua)'
    'Mean_Norm Gait Speed (m.s^{-1})'

    % --- Variabilité ---
    'CV_Single support time (%)'
    'CV_Double support time (%)'
    'CV_BaseOfSupport (cm)'
    'CV_StepWidth (cm)'
    'CV_Gait speed (m.s^{-1})'
    'CV_Stride length (m)'
    'CV_Stride time (s)'
    'CV_WalkRatio'
    'CV_Norm WR (ua)'
    'CV_Cadence (step.min^{-1})'
    'CV_Norm Step length (ua)'
    'CV_Norm Cadence (ua)'
    'CV_Norm StepWidth (ua)'
    'CV_Norm Gait Speed (m.s^{-1})'

    % --- Stabilité dynamique ---
    'Mean_MoS AP HS (mm)'
    'Mean_MoS ML HS (mm)'
    'Mean_MoS AP Stance (mm)'
    'Mean_MoS ML Stance (mm)'
    'Mean_MoS AP HS (%L0)'
    'Mean_MoS ML HS (%L0)'
    'Mean_MoS AP Stance (%L0)'
    'Mean_MoS ML Stance (%L0)'

    % --- Smoothness ---
    'Mean_COM SPARC Magnitude (ua)'
    'Mean_COM LDLJ Magnitude (ua)'
    'Mean_STERN SPARC Magnitude (ua)'
    'Mean_STERN LDLJ Magnitude (ua)'

    % --- Indices de symétrie ---
    'SI_Stride time (s)'
    'SI_Stride length (m)'
    'SI_Double support time (%)'
    'SI_Single support time (%)'
    'SI_BaseOfSupport (cm)'
    'SI_StepWidth (cm)'
    'SI_WalkRatio'
    'SI_Norm WR (ua)'
    'SI_Cadence (step.min^{-1})'
    'SI_Norm Step length (ua)'
    'SI_Norm Cadence (ua)'
    'SI_Norm StepWidth (ua)'
};

%% =============== Préparation des données (format long) ===============
fprintf('=== PRÉPARATION DES DONNÉES & EXPORT DATA FORMAT LONG ===\n');

% Fusion de toutes les conditions
DATA_all = [SpatioTemporalDATA.ALL.Plat;
            SpatioTemporalDATA.ALL.Medium;
            SpatioTemporalDATA.ALL.High];

% Harmoniser le facteur Surface
if ismember('Condition', DATA_all.Properties.VariableNames) && ~ismember('Surface', DATA_all.Properties.VariableNames)
    DATA_all.Properties.VariableNames{'Condition'} = 'Surface';
end

% Sécuriser colonnes minimales
if ~ismember('Participant', DATA_all.Properties.VariableNames)
    error('La colonne "Participant" est introuvable dans DATA_all.');
end
if ~ismember('Surface', DATA_all.Properties.VariableNames)
    error('La colonne "Surface" (ou "Condition") est introuvable dans DATA_all.');
end

% Convertir en categoricals
DATA_all.Participant = categorical(DATA_all.Participant);
DATA_all.Surface     = categorical(DATA_all.Surface, surfaces); % force l'ordre

% Mapping Participant -> AgeGroup
fprintf('Construction du mapping Participant -> AgeGroup...\n');
p2g = containers.Map('KeyType','char','ValueType','char');

for g = 1:numel(groups)
    gName = groups{g};
    if ~isfield(SpatioTemporalDATA, gName), continue; end

    for s = 1:numel(surfaces)
        surfName = surfaces{s};
        if ~isfield(SpatioTemporalDATA.(gName), surfName), continue; end

        T = SpatioTemporalDATA.(gName).(surfName);
        if isempty(T) || ~istable(T) || ~ismember('Participant', T.Properties.VariableNames), continue; end

        parts = unique(string(T.Participant));
        for k = 1:numel(parts)
            key = char(parts(k));
            if ~isKey(p2g, key)
                p2g(key) = gName;
            end
        end
    end
end

% Ajouter AgeGroup
DATA_all.AgeGroup = repmat({''}, height(DATA_all), 1);
parts_all = string(DATA_all.Participant);

for i = 1:numel(parts_all)
    key = char(parts_all(i));
    if isKey(p2g, key)
        DATA_all.AgeGroup{i} = p2g(key);
    else
        DATA_all.AgeGroup{i} = '';
    end
end

% Retirer lignes sans groupe
missingGroup = cellfun(@isempty, DATA_all.AgeGroup);
if any(missingGroup)
    warning('%.0f lignes sans AgeGroup (participants non mappés). Elles seront retirées.', sum(missingGroup));
    DATA_all(missingGroup, :) = [];
end

% Catégorielles avec ordre
DATA_all.AgeGroup = categorical(DATA_all.AgeGroup, groups, 'Ordinal', true);

% Vérifier distribution (Plat)
fprintf('\nDistribution des participants (surface Plat):\n');
for g = 1:numel(groups)
    n = sum(DATA_all.AgeGroup == groups{g} & DATA_all.Surface == 'Plat');
    fprintf('  %s: %d participants\n', groups{g}, n);
end

% =============== Sauvegarde des données préparées (DATA_all) ===============
fprintf('\nSauvegarde de DATA_all...\n');

prep_path = fullfile(stats_path, 'Prepared_Data');
if ~exist(prep_path, 'dir'), mkdir(prep_path); end

save(fullfile(prep_path, 'DATA_all_prepared.mat'), 'DATA_all', '-v7.3');
fprintf('✓ DATA_all_prepared.mat sauvegardé\n');

DATA_all_csv = DATA_all;
varsCat = varfun(@iscategorical, DATA_all_csv, 'OutputFormat', 'uniform');
catNames = DATA_all_csv.Properties.VariableNames(varsCat);
for k = 1:numel(catNames)
    DATA_all_csv.(catNames{k}) = string(DATA_all_csv.(catNames{k}));
end
writetable(DATA_all_csv, fullfile(prep_path, 'DATA_all_prepared.csv'));
fprintf('✓ DATA_all_prepared.csv sauvegardé\n');

%% =============== Analyses LMM pour chaque variable ===============
results_summary  = table();
detailed_results = struct();
failed_analyses  = {};

fprintf('\n=== ANALYSE STATISTIQUE LMM ===\n');

allVarNames = DATA_all.Properties.VariableNames;

for v = 1:numel(variables_to_test)
    label = variables_to_test{v};
    fprintf('\n--- Variable %d/%d: %s ---\n', v, numel(variables_to_test), label);
    yName = resolve_varname(label, allVarNames);

    if isempty(yName)
        warning('Variable non trouvée dans la table (label: "%s").', label);
        failed_analyses{end+1} = label; 
        continue;
    end

    try
        % Préparer les données (enlever NaN)
        cols = {'Participant','AgeGroup','Surface', yName};
        data_clean = DATA_all(:, cols);
        data_clean = rmmissing(data_clean);

        if height(data_clean) < 30
            warning('Pas assez de données pour: %s (n=%d)', label, height(data_clean));
            failed_analyses{end+1} = label; 
            continue;
        end

        % Ajuster le modèle LMM
        [lmm_results, lme, model_meta] = fit_lmm_model(data_clean, yName, FIT_METHOD, DF_METHOD);

        % Diagnostics
        assumptions = check_lmm_assumptions(lme);

        % Stocker les résultats détaillés
        safe_field = matlab.lang.makeValidName(yName);
        detailed_results.(safe_field) = struct( ...
            'Label',       label, ...
            'YName',       yName, ...
            'Model',       lme, ...
            'ModelMeta',   model_meta, ...
            'Results',     lmm_results, ...
            'Assumptions', assumptions);

        % Ajouter au résumé
        new_row = table();
        new_row.Label   = {label};
        new_row.YName   = {yName};
        new_row.N_obs   = height(data_clean);

        new_row.Surface_F   = lmm_results.Surface.F;
        new_row.Surface_DF1 = lmm_results.Surface.DF1;
        new_row.Surface_DF2 = lmm_results.Surface.DF2;
        new_row.Surface_p   = lmm_results.Surface.p;
        new_row.Surface_sig = {get_sig_stars(lmm_results.Surface.p)};

        new_row.AgeGroup_F   = lmm_results.AgeGroup.F;
        new_row.AgeGroup_DF1 = lmm_results.AgeGroup.DF1;
        new_row.AgeGroup_DF2 = lmm_results.AgeGroup.DF2;
        new_row.AgeGroup_p   = lmm_results.AgeGroup.p;
        new_row.AgeGroup_sig = {get_sig_stars(lmm_results.AgeGroup.p)};

        new_row.Interaction_F   = lmm_results.Interaction.F;
        new_row.Interaction_DF1 = lmm_results.Interaction.DF1;
        new_row.Interaction_DF2 = lmm_results.Interaction.DF2;
        new_row.Interaction_p   = lmm_results.Interaction.p;
        new_row.Interaction_sig = {get_sig_stars(lmm_results.Interaction.p)};

        new_row.AIC = lmm_results.AIC;
        new_row.BIC = lmm_results.BIC;
        new_row.LogLik = lmm_results.LogLikelihood;

        % Métadonnées LMM
        new_row.FitMethod = {model_meta.FitMethod};
        new_row.DFMethod  = {model_meta.DFMethod};

        results_summary = [results_summary; new_row]; %#ok<AGROW>

    catch ME
        warning('Erreur pour "%s" (Y="%s"): %s', label, yName, ME.message);
        failed_analyses{end+1} = label;
    end
end

%% =============== Correction comparaisons multiples (FDR) ===============
fprintf('\n=== CORRECTION POUR COMPARAISONS MULTIPLES (FDR) ===\n');

% Architecture FDR:
% - FDR appliquée séparément pour chaque famille d'effets, à travers les variables:
%   (i) Surface p-values, (ii) AgeGroup p-values, (iii) Interaction p-values
if ~isempty(results_summary)
    results_summary = apply_fdr_correction(results_summary);

    n_sig_surface = sum(results_summary.Surface_p_FDR < 0.05, 'omitnan');
    n_sig_age     = sum(results_summary.AgeGroup_p_FDR < 0.05, 'omitnan');
    n_sig_inter   = sum(results_summary.Interaction_p_FDR < 0.05, 'omitnan');

    fprintf('Effets significatifs après correction FDR (par famille d''effets):\n');
    fprintf('  Surface: %d/%d (%.1f%%)\n', n_sig_surface, height(results_summary), 100*n_sig_surface/height(results_summary));
    fprintf('  Groupe d''âge: %d/%d (%.1f%%)\n', n_sig_age, height(results_summary), 100*n_sig_age/height(results_summary));
    fprintf('  Interaction: %d/%d (%.1f%%)\n', n_sig_inter, height(results_summary), 100*n_sig_inter/height(results_summary));

for i = 1:height(results_summary)

    % récupérer la variable
    yName = results_summary.YName{i};
    label = results_summary.Label{i};

    % reconstruire les données
    data_plot = DATA_all(:, {'Participant','AgeGroup','Surface', yName});
    data_plot = rmmissing(data_plot);

    % construire une structure p-values pour la figure
    pvals_fig = struct();
    pvals_fig.Surface_p      = results_summary.Surface_p(i);
    pvals_fig.Surface_p_FDR  = results_summary.Surface_p_FDR(i);
    pvals_fig.AgeGroup_p     = results_summary.AgeGroup_p(i);
    pvals_fig.AgeGroup_p_FDR = results_summary.AgeGroup_p_FDR(i);
    pvals_fig.Inter_p        = results_summary.Interaction_p(i);
    pvals_fig.Inter_p_FDR    = results_summary.Interaction_p_FDR(i);

    plot_lmm_results(data_plot, yName, label, pvals_fig, stats_path);
end

else
    warning('Aucun résultat dans results_summary (tout a échoué ou aucune variable trouvée).');
end

%% =============== Exports ===============
fprintf('\n=== EXPORT DES RÉSULTATS ===\n');

writetable(results_summary, fullfile(stats_path, 'LMM_Summary.csv'));
fprintf('✓ Tableau récapitulatif exporté (CSV)\n');

save(fullfile(stats_path, 'Detailed_LMM_Results.mat'), 'detailed_results');
fprintf('✓ Résultats détaillés sauvegardés (.mat)\n');

generate_lmm_report(results_summary, detailed_results, failed_analyses, stats_path, surfaces, groups);
fprintf('✓ Rapport texte généré\n');

export_to_excel(results_summary, stats_path);
fprintf('✓ Rapport Excel exporté\n');

fprintf('\nANALYSES TERMINÉES.\n');
fprintf('Résultats dans: %s\n', stats_path);
fprintf('Variables analysées: %d/%d\n', height(results_summary), numel(variables_to_test));
if ~isempty(failed_analyses)
    fprintf('Analyses échouées: %d\n', numel(failed_analyses));
end

%% ========================= FONCTIONS LOCALES =========================

function yName = resolve_varname(label, tableVarNames)
    yName = '';

    if any(strcmp(tableVarNames, label))
        yName = label; return;
    end

    cand = matlab.lang.makeValidName(label);
    if any(strcmp(tableVarNames, cand))
        yName = cand; return;
    end

    clean = @(s) lower(regexprep(s, '[^a-zA-Z0-9]+', ''));
    label_c = clean(label);

    tv = tableVarNames;
    tv_clean = cellfun(clean, tv, 'UniformOutput', false);

    idx = find(strcmp(tv_clean, label_c), 1);
    if ~isempty(idx)
        yName = tv{idx}; return;
    end

    idx = find(contains(tv_clean, label_c), 1);
    if ~isempty(idx)
        yName = tv{idx}; return;
    end
end

function [lmm_results, lme, meta] = fit_lmm_model(data, yName, fitMethod, dfMethod)
% - FitMethod: 'REML' ou 'ML'
% - DFMethod:  'Satterthwaite' (si supporté) sinon fallback 'Residual'
% - Interaction testée via coefTest -> pas de NaN

    if nargin < 3 || isempty(fitMethod), fitMethod = 'REML'; end
    if nargin < 4 || isempty(dfMethod),  dfMethod  = 'Satterthwaite'; end

    data.Y = data.(yName);

    formula = 'Y ~ Surface * AgeGroup + (1|Participant)';
    lme = fitlme(data, formula, 'FitMethod', fitMethod);

    % DFMethod: fallback si non supporté
    try
        anova_tbl = anova(lme, 'DFMethod', dfMethod);
        dfMethodUsed = dfMethod;
    catch
        anova_tbl = anova(lme, 'DFMethod', 'Residual');
        dfMethodUsed = 'Residual';
        warning('DFMethod "%s" non supporté -> utilisation de "Residual".', dfMethod);
    end

    meta = struct('Formula', formula, 'FitMethod', fitMethod, 'DFMethod', dfMethodUsed);

    lmm_results = struct();

    function out = getTerm(termName)
        idx = find(strcmp(anova_tbl.Term, termName), 1);
        if isempty(idx)
            out = struct('F', NaN, 'p', NaN, 'DF1', NaN, 'DF2', NaN);
        else
            out = struct( ...
                'F',   anova_tbl.FStat(idx), ...
                'p',   anova_tbl.pValue(idx), ...
                'DF1', anova_tbl.DF1(idx), ...
                'DF2', anova_tbl.DF2(idx));
        end
    end

    % Effets principaux
    lmm_results.Surface  = getTerm('Surface');
    lmm_results.AgeGroup = getTerm('AgeGroup');

    % Interaction: test global via coefTest (Wald F-test)
    betaNames = string(lme.CoefficientNames);
    idxInter = contains(betaNames, 'Surface') & contains(betaNames, ':') & contains(betaNames, 'AgeGroup');

    if any(idxInter)
        cols = find(idxInter);
        L = zeros(numel(cols), numel(betaNames));
        for k = 1:numel(cols)
            L(k, cols(k)) = 1;
        end

        [pInt, FInt, df1Int, df2Int] = coefTest(lme, L);
        lmm_results.Interaction = struct('F', FInt, 'p', pInt, 'DF1', df1Int, 'DF2', df2Int);
    else
        lmm_results.Interaction = struct('F', NaN, 'p', NaN, 'DF1', NaN, 'DF2', NaN);
        warning('Aucun coefficient d''interaction Surface:AgeGroup détecté dans lme.CoefficientNames.');
    end

    % Critères d'information
    lmm_results.AIC = lme.ModelCriterion.AIC;
    lmm_results.BIC = lme.ModelCriterion.BIC;
    lmm_results.LogLikelihood = lme.LogLikelihood;

    % Affichage
    fprintf('  Surface: F(%g,%g)=%.3f, p=%.4f %s\n', ...
        lmm_results.Surface.DF1, lmm_results.Surface.DF2, ...
        lmm_results.Surface.F, lmm_results.Surface.p, get_sig_stars(lmm_results.Surface.p));

    fprintf('  AgeGroup: F(%g,%g)=%.3f, p=%.4f %s\n', ...
        lmm_results.AgeGroup.DF1, lmm_results.AgeGroup.DF2, ...
        lmm_results.AgeGroup.F, lmm_results.AgeGroup.p, get_sig_stars(lmm_results.AgeGroup.p));

    fprintf('  Interaction: F(%g,%g)=%.3f, p=%.4f %s\n', ...
        lmm_results.Interaction.DF1, lmm_results.Interaction.DF2, ...
        lmm_results.Interaction.F, lmm_results.Interaction.p, get_sig_stars(lmm_results.Interaction.p));
end

function assumptions = check_lmm_assumptions(lme)
    assumptions = struct();

    residuals = lme.Residuals.Raw;

    [h_ks, p_ks] = kstest(zscore(residuals));
    assumptions.normality.test = 'Kolmogorov-Smirnov';
    assumptions.normality.p    = p_ks;
    assumptions.normality.ok   = (h_ks == 0);

    fitted = lme.Fitted;
    [~, p_h] = corr(abs(residuals), fitted, 'type', 'Spearman', 'rows', 'complete');
    assumptions.homoscedasticity.test = 'abs(resid) vs fitted (Spearman)';
    assumptions.homoscedasticity.p    = p_h;
    assumptions.homoscedasticity.ok   = (p_h > 0.05);

    std_resid = residuals ./ std(residuals);
    n_out = sum(abs(std_resid) > 3);
    assumptions.outliers.n       = n_out;
    assumptions.outliers.percent = 100 * n_out / numel(residuals);
    assumptions.outliers.ok      = (assumptions.outliers.percent < 5);

    fprintf('  Assumptions:\n');
    fprintf('    Normalité: p=%.4f %s\n', p_ks, okfail(assumptions.normality.ok));
    fprintf('    Homoscédasticité: p=%.4f %s\n', p_h, okfail(assumptions.homoscedasticity.ok));
    fprintf('    Outliers: %.1f%% %s\n', assumptions.outliers.percent, okfail(assumptions.outliers.ok));
end

function s = okfail(tf)
    if tf, s = 'OK'; else, s = 'FAIL'; end
end

function plot_lmm_results(data, yName, label, pvals, save_path)
    surfCats  = categories(data.Surface);
    groupCats = categories(data.AgeGroup);

    fig = figure('Position', [100, 100, 1400, 500]);

    subplot(1,3,1);
    means = zeros(numel(surfCats), 1);
    sems  = zeros(numel(surfCats), 1);
    for i = 1:numel(surfCats)
        idx = (data.Surface == surfCats{i});
        vals = data{idx, yName};
        means(i) = mean(vals, 'omitnan');
        sems(i)  = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
    end
    bar(means); hold on;
    errorbar(1:numel(surfCats), means, sems, 'k.', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', surfCats);
    ylabel(label, 'Interpreter','none');
    title(sprintf('Surface: p = %.4f | p_{FDR} = %.4f', pvals.Surface_p, pvals.Surface_p_FDR), 'Interpreter','tex');
    grid on;

    subplot(1,3,2);
    means = zeros(numel(groupCats), 1);
    sems  = zeros(numel(groupCats), 1);
    for i = 1:numel(groupCats)
        idx = (data.AgeGroup == groupCats{i});
        vals = data{idx, yName};
        means(i) = mean(vals, 'omitnan');
        sems(i)  = std(vals, 'omitnan') / sqrt(sum(~isnan(vals)));
    end
    bar(means); hold on;
    errorbar(1:numel(groupCats), means, sems, 'k.', 'LineWidth', 1.5);
    set(gca, 'XTickLabel', groupCats, 'XTickLabelRotation', 45);
    ylabel(label, 'Interpreter','none');
    title(sprintf('AgeGroup: p = %.4f | p_{FDR} = %.4f', pvals.AgeGroup_p, pvals.AgeGroup_p_FDR), 'Interpreter','tex');
    grid on;

    subplot(1,3,3);
    for g = 1:numel(groupCats)
        mg = zeros(numel(surfCats), 1);
        for s = 1:numel(surfCats)
            idx = (data.Surface == surfCats{s}) & (data.AgeGroup == groupCats{g});
            mg(s) = mean(data{idx, yName}, 'omitnan');
        end
        plot(1:numel(surfCats), mg, '-o', 'LineWidth', 2, 'DisplayName', char(groupCats{g}));
        hold on;
    end
    set(gca, 'XTick', 1:numel(surfCats), 'XTickLabel', surfCats);
    ylabel(label, 'Interpreter','none');
    title(sprintf('Interaction: p = %.4f | p_{FDR} = %.4f', pvals.Inter_p, pvals.Inter_p_FDR), 'Interpreter','tex');
    legend('Location','best'); grid on;

    safe_name = matlab.lang.makeValidName(yName);
    saveas(fig, fullfile(save_path, ['Plot_' safe_name '.png']));
    close(fig);
end

function results_summary = apply_fdr_correction(results_summary)
% FDR Benjamini-Hochberg, appliqué séparément par famille d'effets:
% - Surface_p (toutes les variables)
% - AgeGroup_p (toutes les variables)
% - Interaction_p (toutes les variables)

    [~, ~, ~, pS] = fdr_bh(results_summary.Surface_p);
    results_summary.Surface_p_FDR = pS;
    results_summary.Surface_sig_FDR = cellfun(@get_sig_stars, num2cell(pS), 'UniformOutput', false);

    [~, ~, ~, pG] = fdr_bh(results_summary.AgeGroup_p);
    results_summary.AgeGroup_p_FDR = pG;
    results_summary.AgeGroup_sig_FDR = cellfun(@get_sig_stars, num2cell(pG), 'UniformOutput', false);

    [~, ~, ~, pI] = fdr_bh(results_summary.Interaction_p);
    results_summary.Interaction_p_FDR = pI;
    results_summary.Interaction_sig_FDR = cellfun(@get_sig_stars, num2cell(pI), 'UniformOutput', false);
end

function [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, ~, ~)
    if nargin < 2, q = 0.05; end

    p = pvals(:);
    m = numel(p);

    [ps, idx] = sort(p);
    thresh = (1:m)' * q / m;

    below = ps <= thresh;
    if any(below)
        k = find(below, 1, 'last');
        crit_p = ps(k);
        h = p <= crit_p;
    else
        crit_p = 0;
        h = false(m,1);
    end

    % adjusted p-values (BH)
    adj_sorted = ps .* m ./ (1:m)';
    adj_sorted = flipud(cummin(flipud(adj_sorted)));
    adj_sorted(adj_sorted > 1) = 1;

    adj_p = nan(m,1);
    adj_p(idx) = adj_sorted;

    adj_ci_cvrg = NaN;
end

function stars = get_sig_stars(p)
    if isnan(p)
        stars = 'na';
    elseif p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'ns';
    end
end

function generate_lmm_report(results_summary, detailed_results, failed_analyses, stats_path, surfaces, groups)
    fid = fopen(fullfile(stats_path, 'LMM_Statistical_Report.txt'), 'w');
    if fid < 0
        warning('Impossible de créer le rapport texte.');
        return;
    end

    fprintf(fid, "========================================================\n");
    fprintf(fid, "RAPPORT STATISTIQUE - MODÈLES LINÉAIRES MIXTES (LMM)\n");
    fprintf(fid, "========================================================\n\n");
    fprintf(fid, "Date: %s\n\n", datestr(now));

    fprintf(fid, "DESIGN:\n");
    fprintf(fid, "  Surface (répétée): %s\n", strjoin(surfaces, ', '));
    fprintf(fid, "  Groupe d'âge (inter-sujets): %s\n\n", strjoin(groups, ', '));

    fprintf(fid, "MODÈLE (fixé):\n");
    fprintf(fid, "  Y ~ Surface * AgeGroup + (1|Participant)\n\n");

    % Méta (FitMethod / DFMethod) si dispo
    fns = fieldnames(detailed_results);
    if ~isempty(fns) && isfield(detailed_results.(fns{1}), 'ModelMeta')
        MM = detailed_results.(fns{1}).ModelMeta;
        fprintf(fid, "ESTIMATION / DF:\n");
        fprintf(fid, "  FitMethod: %s\n", MM.FitMethod);
        fprintf(fid, "  DFMethod:  %s\n\n", MM.DFMethod);
    end

    fprintf(fid, "MULTIPLE TESTING:\n");
    fprintf(fid, "  Benjamini–Hochberg FDR applied across outcomes separately for each effect family:\n");
    fprintf(fid, "   - Surface p-values (one per outcome)\n");
    fprintf(fid, "   - AgeGroup p-values (one per outcome)\n");
    fprintf(fid, "   - Surface×AgeGroup interaction p-values (one per outcome)\n\n");

    fprintf(fid, "SEUIL:\n");
    fprintf(fid, "  alpha = 0.05\n\n");

    fprintf(fid, "========================================================\n");
    fprintf(fid, "RÉSULTATS (résumé)\n");
    fprintf(fid, "========================================================\n\n");

    hasFDR = ismember('Surface_p_FDR', results_summary.Properties.VariableNames);

    for i = 1:height(results_summary)
        fprintf(fid, "%d) %s\n", i, results_summary.Label{i});
        fprintf(fid, "   YName: %s\n", results_summary.YName{i});
        fprintf(fid, "   N: %d\n", results_summary.N_obs(i));

        fprintf(fid, "   Surface:     F(%g,%g)=%.3f, p=%.4f %s\n", ...
            results_summary.Surface_DF1(i), results_summary.Surface_DF2(i), ...
            results_summary.Surface_F(i), results_summary.Surface_p(i), ...
            results_summary.Surface_sig{i});

        fprintf(fid, "   AgeGroup:    F(%g,%g)=%.3f, p=%.4f %s\n", ...
            results_summary.AgeGroup_DF1(i), results_summary.AgeGroup_DF2(i), ...
            results_summary.AgeGroup_F(i), results_summary.AgeGroup_p(i), ...
            results_summary.AgeGroup_sig{i});

        fprintf(fid, "   Interaction: F(%g,%g)=%.3f, p=%.4f %s\n", ...
            results_summary.Interaction_DF1(i), results_summary.Interaction_DF2(i), ...
            results_summary.Interaction_F(i), results_summary.Interaction_p(i), ...
            results_summary.Interaction_sig{i});

        if hasFDR
            fprintf(fid, "   FDR: Surface p=%.4f (%s), AgeGroup p=%.4f (%s), Interaction p=%.4f (%s)\n", ...
                results_summary.Surface_p_FDR(i), results_summary.Surface_sig_FDR{i}, ...
                results_summary.AgeGroup_p_FDR(i), results_summary.AgeGroup_sig_FDR{i}, ...
                results_summary.Interaction_p_FDR(i), results_summary.Interaction_sig_FDR{i});
        end

        fprintf(fid, "   Fit: AIC=%.2f, BIC=%.2f, LogLik=%.2f\n\n", ...
            results_summary.AIC(i), results_summary.BIC(i), results_summary.LogLik(i));
    end

    if ~isempty(failed_analyses)
        fprintf(fid, "========================================================\n");
        fprintf(fid, "VARIABLES ÉCHOUÉES / NON TROUVÉES\n");
        fprintf(fid, "========================================================\n");
        for k = 1:numel(failed_analyses)
            fprintf(fid, "- %s\n", failed_analyses{k});
        end
        fprintf(fid, "\n");
    end

    % Diagnostics
    fprintf(fid, "========================================================\n");
    fprintf(fid, "DIAGNOSTICS (aperçu)\n");
    fprintf(fid, "========================================================\n");
    for k = 1:numel(fns)
        S = detailed_results.(fns{k});
        if isfield(S, 'Assumptions')
            A = S.Assumptions;
            fprintf(fid, "- %s:\n", S.Label);
            fprintf(fid, "   Normalité (KS) p=%.4f (%s)\n", A.normality.p, okfail(A.normality.ok));
            fprintf(fid, "   Homoscédasticité p=%.4f (%s)\n", A.homoscedasticity.p, okfail(A.homoscedasticity.ok));
            fprintf(fid, "   Outliers: %.1f%% (%s)\n", A.outliers.percent, okfail(A.outliers.ok));
        end
    end

    fclose(fid);
end

function export_to_excel(results_summary, stats_path)
    excel_file = fullfile(stats_path, 'LMM_Results_Complete.xlsx');

    writetable(results_summary, excel_file, 'Sheet', 'Summary');

    sig_mask = results_summary.Surface_p < 0.05 | results_summary.AgeGroup_p < 0.05 | results_summary.Interaction_p < 0.05;
    if any(sig_mask)
        writetable(results_summary(sig_mask,:), excel_file, 'Sheet', 'Significant_Only');
    else
        writetable(results_summary(1:min(1,height(results_summary)),:), excel_file, 'Sheet', 'Significant_Only');
    end

    model_cols = intersect({'Label','YName','AIC','BIC','LogLik','FitMethod','DFMethod'}, results_summary.Properties.VariableNames, 'stable');
    writetable(results_summary(:, model_cols), excel_file, 'Sheet', 'Model_Fit');
end