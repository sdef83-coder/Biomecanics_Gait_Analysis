%% ANOVA À MESURES RÉPÉTÉES À DEUX FACTEURS
% Facteur 1: Groupe d'âge (between-subject)
% Facteur 2: Surface (within-subject - mesures répétées)
clc, clear, close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\gaitAnalysisGUI\result');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\gaitAnalysisGUI\functions'));
save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\gaitAnalysisGUI\result'

% Charger les données après avoir exécuté votre script principal
load(fullfile(save_path, 'SpatioTemporalDATA.mat'));

%%

% Définir les variables à analyser
variablesToAnalyze = {'Mean_Single support time (%)','Mean_Double support time (%)','Mean_Stride width (mm)', ...
                      'Mean_Gait speed (m.s^{-1})','Mean_Stride length (m)', ...
                      'Mean_Normalized Walk ratio (ua)','Mean_Cadence (step.min^{-1})', ...
                      'Mean_Normalized Step length (ua)', 'Mean_Normalized Cadence (ua)', ...
                      'CV_Single support time (%)','CV_Double support time (%)','CV_Stride width (mm)', ...
                      'CV_Gait speed (m.s^{-1})','CV_Stride length (m)', ...
                      'CV_Normalized Walk ratio (ua)','CV_Cadence (step.min^{-1})', ...
                      'CV_Normalized Step length (ua)', 'CV_Normalized Cadence (ua)'}; % Mean + CV pour 9 Variables

% Groupes d'âge disponibles
groupList = {'JeunesEnfants', 'Enfants', 'Adolescents', 'Adultes'};
Condition = {'Plat', 'Medium', 'High'};

% Résultats ANOVA
anovaResults = struct();

% Fonction utilitaire pour échapper les underscores dans les noms de variables (pour affichage sans erreur)
cleanVarName = @(v) strrep(v, '_', '\_');

% BOUCLE PRINCIPALE POUR CHAQUE VARIABLE
for iVar = 1:length(variablesToAnalyze)
    varName = variablesToAnalyze{iVar};
    fprintf('\n=== ANALYSE ANOVA POUR : %s ===\n', varName);
    
    % Préparer les données pour l'ANOVA
    allData = [];
    participantIDs = {};
    groupLabels = {};
    surfaceLabels = {};
    
    participantCounter = 1;
    
    % Parcourir chaque groupe d'âge
    for g = 1:length(groupList)
        groupName = groupList{g};
        
        if isfield(SpatioTemporalDATA, groupName)
            % Obtenir la liste des participants pour ce groupe
            if isfield(SpatioTemporalDATA.(groupName), 'Plat')
                groupTable = SpatioTemporalDATA.(groupName).Plat;
                if ~isempty(groupTable) && any(strcmp(groupTable.Properties.VariableNames, varName))
                    participants = unique(groupTable.Participant);
                    
                    % Pour chaque participant
                    for p = 1:length(participants)
                        participantID = participants{p};
                        participantData = [];
                        validData = true;
                        
                        % Collecter les données pour les 3 surfaces
                        for iC = 1:length(Condition)
                            cond = Condition{iC};
                            
                            if isfield(SpatioTemporalDATA.(groupName), cond)
                                condTable = SpatioTemporalDATA.(groupName).(cond);
                                
                                if ~isempty(condTable) && any(strcmp(condTable.Properties.VariableNames, varName))
                                    participantRow = strcmp(condTable.Participant, participantID);
                                    
                                    if any(participantRow)
                                        value = condTable.(varName)(participantRow);
                                        if ~isempty(value) && ~isnan(value(1))
                                            participantData(end+1) = value(1);
                                        else
                                            validData = false;
                                            break;
                                        end
                                    else
                                        validData = false;
                                        break;
                                    end
                                else
                                    validData = false;
                                    break;
                                end
                            else
                                validData = false;
                                break;
                            end
                        end
                        
                        % Si le participant a des données complètes pour les 3 surfaces
                        if validData && length(participantData) == 3
                            allData = [allData; participantData];
                            participantIDs{participantCounter} = participantID;
                            groupLabels{participantCounter} = groupName;
                            participantCounter = participantCounter + 1;
                        end
                    end
                end
            end
        end
    end
    
    
    fprintf('Participants avec données complètes : %d\n', size(allData, 1));
    meanAcross = mean(allData, 2);

    if iscell(groupLabels)
    groupLabelsCat = categorical(groupLabels(:));
else
    groupLabelsCat = categorical(cellstr(groupLabels(:)));
end

    T = table( ...
    categorical(participantIDs(:)), ...
    groupLabelsCat, ...
    allData(:,1), ...
    allData(:,2), ...
    allData(:,3), ...
    'VariableNames', {'Participant', 'Groupe', 'Plat', 'Medium', 'High'});

        % STATISTIQUES DESCRIPTIVES
    fprintf('\n--- STATISTIQUES DESCRIPTIVES (format wide) ---\n');
descriptives = table();

for c = 1:length(Condition)
    condName = Condition{c}; % ex: 'Plat', 'Medium', 'High'
    temp = varfun(@mean, T, 'InputVariables', condName, ...
                  'GroupingVariables', 'Groupe');
    temp_std = varfun(@std, T, 'InputVariables', condName, ...
                      'GroupingVariables', 'Groupe');
    temp_n = varfun(@numel, T, 'InputVariables', condName, ...
                    'GroupingVariables', 'Groupe');

    % Fusionner les colonnes
    temp_all = table();
    temp_all.Groupe = temp.Groupe;
    temp_all.Surface = repmat(categorical({condName}), height(temp), 1);
    temp_all.Mean = temp.(['mean_' condName]);
    temp_all.Std = temp_std.(['std_' condName]);
    temp_all.N = temp_n.(['numel_' condName]);

    % Concaténer
    descriptives = [descriptives; temp_all];
end

disp(descriptives);
anovaResults.(matlab.lang.makeValidName(varName)).descriptives = descriptives;

    % Tester la normalité pour chaque surface
normalityPassed = true;
for col = 1:size(allData, 2)
    [swH, swP] = swtest(allData(:, col), 0.05); % Shapiro-Wilk
    if swH == 1
        normalityPassed = false;
        fprintf('⚠️ Non normal sur surface %d (p=%.4f)\n', col, swP);
    end
end

% Tester l'homogénéité des variances entre groupes pour chaque surface
homogeneityPassed = true;
for col = 1:size(allData, 2)
    if ~iscell(groupLabels)
        groupLabels = cellstr(groupLabels);
    end
    p_levene = vartestn(allData(:, col), groupLabels(:), 'TestType', 'LeveneAbsolute', 'Display', 'off');
    if p_levene < 0.05
        homogeneityPassed = false;
        fprintf('⚠️ Levene non ok pour surface %d (p=%.4f)\n', col, p_levene);
    end
end

% Basculer vers tests non paramétriques si violations
if ~normalityPassed || ~homogeneityPassed
    fprintf('➡️ Passage aux tests NON PARAMÉTRIQUES\n');
    resultNonParam = struct();
    
    % Friedman test pour effet de la surface (intra-sujets)
    p_friedman = friedman(allData, 1, 'off');
    fprintf('Friedman (surface) : p = %.4f\n', p_friedman);

    % === POST-HOC FRIEDMAN SI SIGNIFICATIF ===
if p_friedman < 0.05
    p12 = signrank(allData(:,1), allData(:,2)); % Plat vs Medium
    p13 = signrank(allData(:,1), allData(:,3)); % Plat vs High
    p23 = signrank(allData(:,2), allData(:,3)); % Medium vs High

    fprintf('Post-hoc Friedman (paires Surface) :\n');
    fprintf('  Plat vs Medium : p = %.4f\n', p12);
    fprintf('  Plat vs High   : p = %.4f\n', p13);
    fprintf('  Medium vs High : p = %.4f\n', p23);

    resultNonParam.posthoc_surface = table( ...
        {'Plat vs Medium'; 'Plat vs High'; 'Medium vs High'}, ...
        [p12; p13; p23], ...
        'VariableNames', {'Comparison', 'pValue'});
end

    % Kruskal-Wallis pour effet du groupe sur la moyenne
    meanAcross = mean(allData, 2);
    [p_kruskal, tbl_kw, stats_kw] = kruskalwallis(meanAcross, groupLabels, 'off');
    fprintf('Kruskal-Wallis (groupe) : p = %.4f\n', p_kruskal);

    % === POST-HOC KRUSKAL-WALLIS SI SIGNIFICATIF ===
if p_kruskal < 0.05
    c_groupe_np = multcompare(stats_kw, 'Display', 'off');
    disp('Post-hoc Kruskal-Wallis (groupes) :');
    disp(c_groupe_np);
    resultNonParam.posthoc_groupe = c_groupe_np;
end


    % Sauvegarder dans les résultats
    resultNonParam.p_friedman = p_friedman;
    resultNonParam.p_kruskal = p_kruskal;
    anovaResults.(matlab.lang.makeValidName(varName)).nonParametric = resultNonParam;
    
    continue; % Passer au suivant
end
    
   % PRÉPARATION DES DONNÉES POUR L'ANOVA
    % Reformater les données en format long
    nParticipants = size(allData, 1);
    nSurfaces = size(allData, 2);
    
 % ANOVA À MESURES RÉPÉTÉES (APPROCHE JASP - fitrm + ranova)
try
if iscell(groupLabels)
    groupLabelsCat = categorical(groupLabels(:));
else
    groupLabelsCat = categorical(cellstr(groupLabels(:)));
end

    % Définir le plan intra-sujets (within-subjects)
    within = table(categorical({'Plat'; 'Medium'; 'High'}), 'VariableNames', {'Surface'});

    % Modèle à mesures répétées
    rm = fitrm(T, 'Plat-High ~ Groupe', 'WithinDesign', within);

    % ANOVA à mesures répétées
    ranovatbl = ranova(rm);
    disp(ranovatbl);
    statsSphericity = mauchly(rm); % Tester la sphéricité

    % Calculer la moyenne des 3 surfaces pour chaque participant
T.MeanAcrossSurfaces = mean(T{:, {'Plat', 'Medium', 'High'}}, 2);

% Tester la normalité de la moyenne inter-surface (pour effet du groupe)
[swH_group, swP_group] = swtest(meanAcross, 0.05);
groupNormalityPassed = (swH_group == 0); % true si normalité respectée

if groupNormalityPassed
    % === ANOVA entre-sujets ===
    [p_groupe_tbl, tbl_stats, stats] = anova1(meanAcross, groupLabels, 'off');
    disp(tbl_stats);
    anovaResults.(matlab.lang.makeValidName(varName)).groupMainEffect = tbl_stats;
    anovaResults.(matlab.lang.makeValidName(varName)).groupMainEffect_p = p_groupe_tbl;

    % Post-hoc si significatif
    if p_groupe_tbl < 0.05
        c_groupe = multcompare(stats, 'Display', 'off');
        anovaResults.(matlab.lang.makeValidName(varName)).posthoc_groupe = c_groupe;
        fprintf('Post-hoc Groupe (entre-sujets) :\n'); disp(c_groupe);
    end
else
    % === Test de Kruskal-Wallis ===
    [p_kruskal, tbl_kw, stats_kw] = kruskalwallis(meanAcross, groupLabels, 'off');
    fprintf('⚠️ Test non paramétrique pour Groupe : Kruskal-Wallis, p = %.4f\n', p_kruskal);
    anovaResults.(matlab.lang.makeValidName(varName)).groupMainEffect_p = p_kruskal;
    anovaResults.(matlab.lang.makeValidName(varName)).groupMainEffect = tbl_kw;

    % Post-hoc Kruskal-Wallis si significatif
    if p_kruskal < 0.05
        c_groupe_np = multcompare(stats_kw, 'Display', 'off');
        anovaResults.(matlab.lang.makeValidName(varName)).posthoc_groupe = c_groupe_np;
        fprintf('Post-hoc Groupe (Kruskal-Wallis) :\n'); disp(c_groupe_np);
    end
end

fprintf('Test de normalité (moyenne inter-surfaces) pour l''effet Groupe : p = %.4f -> %s\n', ...
    swP_group, ternary(groupNormalityPassed, 'Normalité OK', '⚠️ Non normal'));

    % Sauvegarder les résultats
    anovaResults.(matlab.lang.makeValidName(varName)).rm = rm;
    anovaResults.(matlab.lang.makeValidName(varName)).ranova = ranovatbl;
    anovaResults.(matlab.lang.makeValidName(varName)).dataWide = T;

    % POST-HOC SI SIGNIFICATIF
    fprintf('\n--- POST-HOC (p < 0.05) ---\n');
    if any(ranovatbl.pValue < 0.05)
        % Effet principal Surface
        if any(strcmp(ranovatbl.Properties.RowNames, '(Intercept):Surface')) && ...
                ranovatbl.pValue(strcmp(ranovatbl.Properties.RowNames, '(Intercept):Surface')) < 0.05
            c_surface = multcompare(rm, 'Surface');
            anovaResults.(matlab.lang.makeValidName(varName)).posthoc_surface = c_surface;
            fprintf('Post-hoc Surface :\n'); disp(c_surface);
        end

        % Interaction Groupe*Surface
        if any(strcmp(ranovatbl.Properties.RowNames, 'Groupe:Surface')) && ...
           ranovatbl.pValue(strcmp(ranovatbl.Properties.RowNames, 'Groupe:Surface')) < 0.05
            c_interaction = multcompare(rm, 'Surface', 'By', 'Groupe');
            anovaResults.(matlab.lang.makeValidName(varName)).posthoc_interaction = c_interaction;
            fprintf('Post-hoc Interaction Groupe × Surface :\n'); disp(c_interaction);
        end

        % Effet principal Groupe
        if isfield(anovaResults.(matlab.lang.makeValidName(varName)), 'groupMainEffect_p') && ...
           anovaResults.(matlab.lang.makeValidName(varName)).groupMainEffect_p < 0.05
            c_groupe = multcompare(stats, 'Display', 'off');
            anovaResults.(matlab.lang.makeValidName(varName)).posthoc_groupe = c_groupe;
            fprintf('Post-hoc Groupe (entre-sujets) :\n'); disp(c_groupe);
        end
    end

catch ME
    fprintf('⚠️ Erreur pour %s : %s\n', varName, ME.message);
    continue;
end

end

% SAUVEGARDER LES RÉSULTATS
save(fullfile(save_path, 'ANOVA_Results.mat'), 'anovaResults');

% RÉSUMÉ FINAL DES RÉSULTATS ANOVA (p < 0.05) + eta²
fprintf('\n=== RÉSUMÉ FINAL DES RÉSULTATS ANOVA (p < 0.05) AVEC η² partiel ou tests non paramétriques ===\n');
varNames = fieldnames(anovaResults);

for i = 1:length(varNames)
    varName = varNames{i};
    fprintf('%s\n', strrep(varName, '_', ' '));

    if isfield(anovaResults.(varName), 'ranova')
        ranovatbl = anovaResults.(varName).ranova;
        statsSphericity = mauchly(anovaResults.(varName).rm);
idxSph = strcmp(statsSphericity.Properties.RowNames, '(Intercept):Surface');

if any(idxSph) && statsSphericity.pValue(idxSph) < 0.05
            fprintf('⚠️ Sphéricité violée (p = %.4f), correction GG appliquée\n', statsSphericity.pValue(idxSph));
        else
            fprintf('✅ Sphéricité respectée (p = %.4f)\n', statsSphericity.pValue(idxSph));
        end

        % Initialisation
        pGroup = NaN; eta2Group = NaN;
        pSurface = NaN; eta2Surface = NaN;
        pInteraction = NaN; eta2Interaction = NaN;

        % Groupe (ANOVA entre-sujets)
        if isfield(anovaResults.(varName), 'groupMainEffect_p') && isfield(anovaResults.(varName), 'groupMainEffect')
            pGroup = anovaResults.(varName).groupMainEffect_p;
            SS_effect = anovaResults.(varName).groupMainEffect{2,2}; % SS groupe
            SS_error  = anovaResults.(varName).groupMainEffect{3,2}; % SS erreur
            eta2Group = SS_effect / (SS_effect + SS_error);
        end

        % Surface (avec correction GG)
        if any(strcmp(ranovatbl.Properties.RowNames, '(Intercept):Surface'))
            idx = strcmp(ranovatbl.Properties.RowNames, '(Intercept):Surface');
            pSurface = ranovatbl.pValueGG(idx); 
            SS_surface = ranovatbl.SumSq(idx);
            SS_error_surface = ranovatbl.SumSq(strcmp(ranovatbl.Properties.RowNames, 'Error(Surface)'));
            eta2Surface = SS_surface / (SS_surface + SS_error_surface);
        end

        % Interaction Groupe×Surface
        if any(strcmp(ranovatbl.Properties.RowNames, 'Groupe:Surface'))
            idx = strcmp(ranovatbl.Properties.RowNames, 'Groupe:Surface');
            pInteraction = ranovatbl.pValueGG(idx); 
            SS_inter = ranovatbl.SumSq(idx);
            SS_error_inter = ranovatbl.SumSq(strcmp(ranovatbl.Properties.RowNames, 'Error(Surface)'));
            eta2Interaction = SS_inter / (SS_inter + SS_error_inter);
        end

        % Affichage avec étoiles
        stars = @(p) repmat('*',1,sum(p < [0.05 0.01 0.001]));
        fprintf('   Groupe   : p = %.4f %s  | η²p = %.3f\n', pGroup, stars(pGroup), eta2Group);
        fprintf('   Surface  : p = %.4f %s  | η²p = %.3f\n', pSurface, stars(pSurface), eta2Surface);
        fprintf('   Interaction Groupe×Surface : p = %.4f %s  | η²p = %.3f\n\n', pInteraction, stars(pInteraction), eta2Interaction);

    elseif isfield(anovaResults.(varName), 'nonParametric')
        np = anovaResults.(varName).nonParametric;
        stars = @(p) repmat('*',1,sum(p < [0.05 0.01 0.001]));
        fprintf('⚠️ Données non normales ou variances inégales : tests NON PARAMÉTRIQUES utilisés\n');
        fprintf('   Kruskal-Wallis (groupe)   : p = %.4f %s\n', np.p_kruskal, stars(np.p_kruskal));
        fprintf('   Friedman (surface)        : p = %.4f %s\n\n', np.p_friedman, stars(np.p_friedman));
    else
        fprintf('⚠️ Aucune analyse disponible pour cette variable.\n\n');
    end
end

% === EXPORT CSV GLOBAL DES RÉSULTATS ANOVA ===
recapTable = table();

for i = 1:length(varNames)
    varName = varNames{i};
    shortName = strrep(varName, '_', ' ');  % Pour un nom lisible

    % Initialisation
    pGroup = NaN; eta2Group = NaN; starsGroup = '';
    pSurface = NaN; eta2Surface = NaN; starsSurface = '';
    pInteraction = NaN; eta2Interaction = NaN; starsInteraction = '';
    testUsed = 'Parametric';

    if isfield(anovaResults.(varName), 'ranova')
        ranovatbl = anovaResults.(varName).ranova;
        
        % Effet Groupe (entre-sujets)
        if isfield(anovaResults.(varName), 'groupMainEffect_p') && isfield(anovaResults.(varName), 'groupMainEffect')
            pGroup = anovaResults.(varName).groupMainEffect_p;
            SS_g = anovaResults.(varName).groupMainEffect{2,2};
            SS_g_err = anovaResults.(varName).groupMainEffect{3,2};
            eta2Group = SS_g / (SS_g + SS_g_err);
            starsGroup = repmat('*', 1, sum(pGroup < [0.05, 0.01, 0.001]));
        end

        % Effet Surface
        if any(strcmp(ranovatbl.Properties.RowNames, '(Intercept):Surface'))
            idx = strcmp(ranovatbl.Properties.RowNames, '(Intercept):Surface');
            pSurface = ranovatbl.pValueGG(idx);
            SS_surf = ranovatbl.SumSq(idx);
            SS_surf_err = ranovatbl.SumSq(strcmp(ranovatbl.Properties.RowNames, 'Error(Surface)'));
            eta2Surface = SS_surf / (SS_surf + SS_surf_err);
            starsSurface = repmat('*', 1, sum(pSurface < [0.05, 0.01, 0.001]));
        end

        % Interaction
        if any(strcmp(ranovatbl.Properties.RowNames, 'Groupe:Surface'))
            idx = strcmp(ranovatbl.Properties.RowNames, 'Groupe:Surface');
            pInteraction = ranovatbl.pValueGG(idx);
            SS_int = ranovatbl.SumSq(idx);
            SS_int_err = ranovatbl.SumSq(strcmp(ranovatbl.Properties.RowNames, 'Error(Surface)'));
            eta2Interaction = SS_int / (SS_int + SS_int_err);
            starsInteraction = repmat('*', 1, sum(pInteraction < [0.05, 0.01, 0.001]));
        end

    elseif isfield(anovaResults.(varName), 'nonParametric')
        testUsed = 'Non-parametric';
        np = anovaResults.(varName).nonParametric;
        pGroup = np.p_kruskal;
        pSurface = np.p_friedman;
        starsGroup = repmat('*', 1, sum(pGroup < [0.05, 0.01, 0.001]));
        starsSurface = repmat('*', 1, sum(pSurface < [0.05, 0.01, 0.001]));
    end

    % Ajouter la ligne au tableau
    recapTable = [recapTable; table({shortName}, {testUsed}, ...
        pGroup, {starsGroup}, eta2Group, ...
        pSurface, {starsSurface}, eta2Surface, ...
        pInteraction, {starsInteraction}, eta2Interaction, ...
        'VariableNames', {'Variable', 'Test', ...
                          'p_Groupe', 'Stars_Groupe', 'Eta2_Groupe', ...
                          'p_Surface', 'Stars_Surface', 'Eta2_Surface', ...
                          'p_Interaction', 'Stars_Interaction', 'Eta2_Interaction'})];
end

% === EXPORTER VERS UN FICHIER CSV ===
outputDir = fullfile(save_path, 'Anova_Recap');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

writetable(recapTable, fullfile(outputDir, 'ANOVA_Recap_Global.csv'));


% RÉSUMÉ FINAL
fprintf('\n=== RÉSUMÉ FINAL ===\n');
fprintf('Analyses ANOVA réalisées pour %d variables\n', length(fieldnames(anovaResults)));
fprintf('Résultats sauvegardés dans : %s\n', fullfile(save_path, 'ANOVA_Results.mat'));
fprintf('\n✅ Fichier CSV global exporté dans : %s\n', fullfile(outputDir, 'ANOVA_Recap_Global.csv'));