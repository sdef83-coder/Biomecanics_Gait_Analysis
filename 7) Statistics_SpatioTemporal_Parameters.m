%% ANOVA À MESURES RÉPÉTÉES À DEUX FACTEURS
% Facteur 1: Groupe d'âge (between-subject)
% Facteur 2: Surface (within-subject - mesures répétées)

clc; clear; close all;

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));

save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Statistics';
if ~exist(save_path, 'dir'); mkdir(save_path); end

% Charger les données après avoir exécuté votre script principal
load('SpatioTemporalDATA.mat');

% Variables issues du main (renommées)
variablesToAnalyze = {
    % --- Moyennes spatio-temporelles ---
    'Mean_Single support time (%)'
    'Mean_Double support time (%)'
    'Mean_Stride width (cm)'
    'Mean_Gait speed (m.s^{-1})'
    'Mean_Stride length (m)'
    'Mean_Stride time (s)'
    'Mean_Norm WR (ua)'
    'Mean_Cadence (step.min^{-1})'
    'Mean_Norm Step length (ua)'
    'Mean_Norm Cadence (ua)'

    % --- Variabilité ---
    'CV_Single support time (%)'
    'CV_Double support time (%)'
    'CV_Stride width (cm)'
    'CV_Gait speed (m.s^{-1})'
    'CV_Stride length (m)'
    'CV_Stride time (s)'
    'CV_Norm WR (ua)'
    'CV_Cadence (step.min^{-1})'
    'CV_Norm Step length (ua)'
    'CV_Norm Cadence (ua)'

    % --- MoS bruts (mm) ---
    'Mean_MoS AP HS (mm)'
    'Mean_MoS ML HS (mm)'
    'Mean_MoS AP Stance (mm)'
    'Mean_MoS ML Stance (mm)'

    % --- MoS normalisés (%L0) ---
    'Mean_MoS AP HS (%L0)'
    'Mean_MoS ML HS (%L0)'
    'Mean_MoS AP Stance (%L0)'
    'Mean_MoS ML Stance (%L0)'

    % --- Indices de symétrie ---
    'SI_Stride time (s)'
    'SI_Stride length (m)'
    'SI_Double support time (%)'
};

groupList = {'JeunesEnfants', 'Enfants', 'Adolescents', 'Adultes'};
Condition = {'Plat', 'Medium', 'High'};

anovaResults = struct();

% BOUCLE PRINCIPALE POUR CHAQUE VARIABLE
for iVar = 1:length(variablesToAnalyze)
    varName = variablesToAnalyze{iVar};
    fprintf('\n=== ANALYSE ANOVA POUR : %s ===\n', varName);
    
    % Préparer les données pour l'ANOVA
    allData = [];
    participantIDs = {};
    groupLabels = {};
    participantCounter = 1;
    
    % Parcourir chaque groupe d'âge
    for g = 1:length(groupList)
        groupName = groupList{g};
        
        if isfield(SpatioTemporalDATA, groupName)
            % Obtenir la liste des participants pour ce groupe à partir de Plat
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
                            participantIDs{participantCounter,1} = participantID;
                            groupLabels{participantCounter,1} = groupName;
                            participantCounter = participantCounter + 1;
                        end
                    end
                end
            end
        end
    end
    
    fprintf('Participants avec données complètes : %d\n', size(allData, 1));
    meanAcross = mean(allData, 2);

    % Table wide
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
        condName = Condition{c};
        temp     = varfun(@mean, T, 'InputVariables', condName, 'GroupingVariables', 'Groupe');
        temp_std = varfun(@std,  T, 'InputVariables', condName, 'GroupingVariables', 'Groupe');
        temp_n   = varfun(@numel,T, 'InputVariables', condName, 'GroupingVariables', 'Groupe');

        temp_all = table();
        temp_all.Groupe  = temp.Groupe;
        temp_all.Surface = repmat(categorical({condName}), height(temp), 1);
        temp_all.Mean    = temp.(['mean_' condName]);
        temp_all.Std     = temp_std.(['std_' condName]);
        temp_all.N       = temp_n.(['numel_' condName]);
        descriptives     = [descriptives; temp_all];
    end
    disp(descriptives);

    % --- Tests d'hypothèses (normalité / homoscédasticité) ---
    normalityPassed    = true;
    homogeneityPassed  = true;

    for col = 1:size(allData, 2)
        [swH, swP] = swtest(allData(:, col), 0.05);
        if swH == 1
            normalityPassed = false;
            fprintf('⚠️ Non normal sur surface %d (p=%.4f)\n', col, swP);
        end

        % Levene
        if ~iscell(groupLabels)
            groupLabels = cellstr(groupLabels);
        end
        p_levene = vartestn(allData(:, col), groupLabels(:), 'TestType', 'LeveneAbsolute', 'Display', 'off');
        if p_levene < 0.05
            homogeneityPassed = false;
            fprintf('⚠️ Levene non ok pour surface %d (p=%.4f)\n', col, p_levene);
        end
    end

    % =========================================================
    % ===============  CAS NON PARAMÉTRIQUE  ==================
    % =========================================================
    if ~normalityPassed || ~homogeneityPassed
        fprintf('➡️ Passage aux tests NON PARAMÉTRIQUES\n');

        % Friedman (surface)
        p_friedman = friedman(allData, 1, 'off');
        fprintf('Friedman (surface) : p = %.4f\n', p_friedman);

        posthoc_surface = [];
        if p_friedman < 0.05
            p12 = signrank(allData(:,1), allData(:,2)); % Plat vs Medium
            p13 = signrank(allData(:,1), allData(:,3)); % Plat vs High
            p23 = signrank(allData(:,2), allData(:,3)); % Medium vs High

            fprintf('Post-hoc Friedman (paires Surface) :\n');
            fprintf('  Plat vs Medium : p = %.4f\n', p12);
            fprintf('  Plat vs High   : p = %.4f\n', p13);
            fprintf('  Medium vs High : p = %.4f\n', p23);

            posthoc_surface = table( ...
                {'Plat vs Medium'; 'Plat vs High'; 'Medium vs High'}, ...
                [p12; p13; p23], ...
                'VariableNames', {'Comparison', 'pValue'});
        end

        % Kruskal-Wallis (groupe)
        meanAcross = mean(allData, 2);
        [p_kruskal, tbl_kw, stats_kw] = kruskalwallis(meanAcross, groupLabels, 'off');
        fprintf('Kruskal-Wallis (groupe) : p = %.4f\n', p_kruskal);

        posthoc_groupe = [];
        if p_kruskal < 0.05
            c_groupe_np = multcompare(stats_kw, 'Display', 'off');
            disp('Post-hoc Kruskal-Wallis (groupes) :');
            disp(c_groupe_np);
            posthoc_groupe = c_groupe_np;
        end

        % ===== stockage propre =====
        resNP = struct();
        resNP.analysisType = 'nonparametric';
        resNP.dataWide     = T;
        resNP.meta.groupsOrder   = groupList;
        resNP.meta.surfacesOrder = Condition;
        resNP.descriptives = descriptives;

        resNP.surface.test  = 'Friedman';
        resNP.surface.p     = p_friedman;
        resNP.surface.posthoc = posthoc_surface;

        resNP.group.test    = 'KruskalWallis';
        resNP.group.p       = p_kruskal;
        resNP.group.posthoc = posthoc_groupe;

        anovaResults.(matlab.lang.makeValidName(varName)) = resNP;
        continue;   % variable suivante
    end

    % =========================================================
    % ===============  CAS PARAMÉTRIQUE  ======================
    % =========================================================
    try
        if iscell(groupLabels)
            groupLabelsCat = categorical(groupLabels(:));
        else
            groupLabelsCat = categorical(cellstr(groupLabels(:)));
        end

        % Design intra-sujet
        within = table(categorical({'Plat'; 'Medium'; 'High'}), 'VariableNames', {'Surface'});

        % fitrm : colonnes de Plat à High
        rm = fitrm(T, 'Plat-High ~ Groupe', 'WithinDesign', within);

        % ANOVA RM
        ranovatbl = ranova(rm);
        disp(ranovatbl);
        statsSphericity = mauchly(rm);

        % moyenne inter-surface pour test groupe
        T.MeanAcrossSurfaces = mean(T{:, {'Plat','Medium','High'}}, 2);
        [swH_group, swP_group] = swtest(meanAcross, 0.05);
        groupNormalityPassed   = (swH_group == 0);

        % --- effet groupe (entre-sujets)
        if groupNormalityPassed
            [p_groupe_tbl, tbl_stats, stats] = anova1(meanAcross, groupLabels, 'off');
        else
            [p_groupe_tbl, tbl_stats, stats] = deal(NaN);
        end

        % ===== stockage propre =====
        res = struct();
        res.analysisType = 'parametric';
        res.dataWide     = T;
        res.meta.groupsOrder   = groupList;
        res.meta.surfacesOrder = Condition;
        res.descriptives = descriptives;
        res.ranova       = ranovatbl;
        res.rm           = rm;
        res.sphericity   = statsSphericity;

        % ----- Surface -----
        res.surface = struct();
        if any(strcmp(ranovatbl.Properties.RowNames, '(Intercept):Surface'))
            idxSurf = strcmp(ranovatbl.Properties.RowNames, '(Intercept):Surface');
            res.surface.test = 'RM-ANOVA (ranova)';
            res.surface.p    = ranovatbl.pValueGG(idxSurf);
            if res.surface.p < 0.05
                res.surface.posthoc = multcompare(rm, 'Surface');
            else
                res.surface.posthoc = [];
            end
        else
            res.surface.test = 'RM-ANOVA (ranova)';
            res.surface.p    = NaN;
            res.surface.posthoc = [];
        end

        % ----- Interaction -----
        res.interaction = struct();
        if any(strcmp(ranovatbl.Properties.RowNames, 'Groupe:Surface'))
            idxInt = strcmp(ranovatbl.Properties.RowNames, 'Groupe:Surface');
            res.interaction.test = 'RM-ANOVA (ranova)';
            res.interaction.p    = ranovatbl.pValueGG(idxInt);
            if res.interaction.p < 0.05
                res.interaction.posthoc = multcompare(rm, 'Surface', 'By','Groupe');
            else
                res.interaction.posthoc = [];
            end
        else
            res.interaction.test = 'RM-ANOVA (ranova)';
            res.interaction.p    = NaN;
            res.interaction.posthoc = [];
        end

        % ----- Groupe (entre-sujets) -----
        res.group = struct();
        if groupNormalityPassed
            res.group.test = 'ANOVA1';
            res.group.p    = p_groupe_tbl;
            res.group.table = tbl_stats;
            if p_groupe_tbl < 0.05
                res.group.posthoc = multcompare(stats, 'Display','off');
            else
                res.group.posthoc = [];
            end
        else
            % fallback non param si tu veux, sinon NaN
            res.group.test = 'ANOVA1 (non normal)';
            res.group.p    = NaN;
            res.group.posthoc = [];
        end

        % Sauvegarde dans struct global
        anovaResults.(matlab.lang.makeValidName(varName)) = res;

    catch ME
        fprintf('⚠️ Erreur pour %s : %s\n', varName, ME.message);
        continue;
    end

end % fin boucle variables

% SAUVEGARDER LES RÉSULTATS
save(fullfile(save_path, 'ANOVA_Results.mat'), 'anovaResults');

% RÉSUMÉ FINAL + EXPORT CSV
fprintf('\n=== RÉSUMÉ FINAL DES RÉSULTATS ANOVA (p < 0.05) AVEC η² partiel OU tests non paramétriques ===\n');

varNames = fieldnames(anovaResults);
recapTable = table();

for i = 1:length(varNames)
    varName = varNames{i};
    fprintf('%s\n', strrep(varName, '_', ' '));

    % init
    pGroup = NaN; eta2Group = NaN; starsGroup = '';
    pSurface = NaN; eta2Surface = NaN; starsSurface = '';
    pInteraction = NaN; eta2Interaction = NaN; starsInteraction = '';
    testUsed = 'Parametric';

    resVar = anovaResults.(varName);

    if isfield(resVar, 'analysisType') && strcmp(resVar.analysisType, 'parametric')
        ranovatbl = resVar.ranova;

        % --- Groupe (entre-sujets)
        if isfield(resVar, 'group') && isfield(resVar.group, 'p')
            pGroup = resVar.group.p;
            starsGroup = repmat('*', 1, sum(pGroup < [0.05 0.01 0.001]));
            % on ne peut pas toujours reconstituer eta2 ici car on a ANOVA1 séparée
        end

        % --- Surface
        if isfield(resVar, 'surface')
            pSurface = resVar.surface.p;
            starsSurface = repmat('*', 1, sum(pSurface < [0.05 0.01 0.001]));
        end

        % --- Interaction
        if isfield(resVar, 'interaction')
            pInteraction = resVar.interaction.p;
            starsInteraction = repmat('*', 1, sum(pInteraction < [0.05 0.01 0.001]));
        end

    elseif isfield(resVar, 'analysisType') && strcmp(resVar.analysisType, 'nonparametric')
        testUsed = 'Non-parametric';
        pGroup   = resVar.group.p;
        pSurface = resVar.surface.p;
        starsGroup   = repmat('*', 1, sum(pGroup < [0.05 0.01 0.001]));
        starsSurface = repmat('*', 1, sum(pSurface < [0.05 0.01 0.001]));
    end

    % ligne recap
    recapTable = [recapTable; table({strrep(varName,'_',' ')}, {testUsed}, ...
        pGroup, {starsGroup}, eta2Group, ...
        pSurface, {starsSurface}, eta2Surface, ...
        pInteraction, {starsInteraction}, eta2Interaction, ...
        'VariableNames', {'Variable', 'Test', ...
                          'p_Groupe', 'Stars_Groupe', 'Eta2_Groupe', ...
                          'p_Surface', 'Stars_Surface', 'Eta2_Surface', ...
                          'p_Interaction', 'Stars_Interaction', 'Eta2_Interaction'})];
end

% Export CSV
outputDir = fullfile(save_path, 'Anova_Recap');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end
writetable(recapTable, fullfile(outputDir, 'ANOVA_Recap_Global.csv'));

fprintf('\n=== RÉSUMÉ FINAL ===\n');
fprintf('Analyses ANOVA réalisées pour %d variables\n', length(fieldnames(anovaResults)));
fprintf('Résultats sauvegardés dans : %s\n', fullfile(save_path, 'ANOVA_Results.mat'));
fprintf('✅ Fichier CSV global exporté dans : %s\n', fullfile(outputDir, 'ANOVA_Recap_Global.csv'));