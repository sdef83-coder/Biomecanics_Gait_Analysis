%% ANOVA À MESURES RÉPÉTÉES À DEUX FACTEURS - VERSION CORRIGÉE
% Facteur 1: Groupe d'âge (between-subject)
% Facteur 2: Surface (within-subject - mesures répétées)

clc; clear; close all;

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));

save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Statistics';
if ~exist(save_path, 'dir'); mkdir(save_path); end

% Charger les données
load('SpatioTemporalDATA.mat');

% Variables à analyser
variablesToAnalyze = {
    % --- Moyennes spatio-temporelles ---
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

groupList = {'JeunesEnfants', 'Enfants', 'Adolescents', 'Adultes'};
groupOrder = categorical(groupList, groupList, 'Ordinal', true);
Condition = {'Plat', 'Medium', 'High'};

anovaResults = struct();

%% BOUCLE PRINCIPALE POUR CHAQUE VARIABLE
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
                                            participantData(end+1) = value(1); %#ok<SAGROW>
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
                            allData = [allData; participantData]; %#ok<AGROW>
                            participantIDs{participantCounter,1} = participantID; %#ok<SAGROW>
                            groupLabels{participantCounter,1} = groupName; %#ok<SAGROW>
                            participantCounter = participantCounter + 1;
                        end
                    end
                end
            end
        end
    end
    
    fprintf('Participants avec données complètes : %d\n', size(allData, 1));

    % *** Garde-fou : si aucune donnée complète, on saute cette variable ***
    if isempty(allData)
        warning('Aucune donnée complète pour la variable %s → variable ignorée.', varName);
        continue;
    end
    
    % *** Vérifier qu'il y a au moins 2 participants par groupe ***
    uniqueGroups = unique(groupLabels);
    groupCounts = zeros(length(uniqueGroups), 1);
    for gg = 1:length(uniqueGroups)
        groupCounts(gg) = sum(strcmp(groupLabels, uniqueGroups{gg}));
    end
    
    if any(groupCounts < 2)
        warning('Au moins un groupe a moins de 2 participants pour %s → variable ignorée.', varName);
        continue;
    end

    meanAcross = mean(allData, 2);

    % Table wide
    groupLabelsCat = categorical(groupLabels(:), groupList, 'Ordinal', true);
    T = table( ...
        categorical(participantIDs(:)), ...
        groupLabelsCat, ...
        allData(:,1), ...
        allData(:,2), ...
        allData(:,3), ...
        'VariableNames', {'Participant', 'Groupe', 'Plat', 'Medium', 'High'});

    % STATISTIQUES DESCRIPTIVES
    fprintf('\n--- STATISTIQUES DESCRIPTIVES ---\n');
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
        descriptives     = [descriptives; temp_all]; %#ok<AGROW>
    end
    disp(descriptives);

    % TESTS D'HYPOTHÈSES
    normalityPassed    = true;
    homogeneityPassed  = true;

    for col = 1:size(allData, 2)
        [swH, swP] = swtest(allData(:, col), 0.05);
        if swH == 1
            normalityPassed = false;
            fprintf('⚠️ Non normalité sur %s (p=%.4f)\n', Condition{col}, swP);
        end

        p_levene = vartestn(allData(:, col), groupLabels(:), ...
            'TestType', 'LeveneAbsolute', 'Display', 'off');
        if p_levene < 0.05
            homogeneityPassed = false;
            fprintf('⚠️ Hétérogénéité des variances sur %s (p=%.4f)\n', Condition{col}, p_levene);
        end
    end

    % =========================================================
    % ===============  CAS NON PARAMÉTRIQUE  ==================
    % =========================================================
    if ~normalityPassed || ~homogeneityPassed
        fprintf('➡️ Analyse NON PARAMÉTRIQUE\n');

        % Friedman (effet Surface - intra-sujet)
        p_friedman = friedman(allData, 1, 'off');
        fprintf('Friedman (Surface) : p = %.4f\n', p_friedman);

        posthoc_surface = [];
        if p_friedman < 0.05
            % Post-hoc avec Wilcoxon signed-rank
            p12 = signrank(allData(:,1), allData(:,2)); % Plat vs Medium
            p13 = signrank(allData(:,1), allData(:,3)); % Plat vs High
            p23 = signrank(allData(:,2), allData(:,3)); % Medium vs High

            fprintf('Post-hoc Friedman :\n');
            fprintf('  Plat vs Medium : p = %.4f\n', p12);
            fprintf('  Plat vs High   : p = %.4f\n', p13);
            fprintf('  Medium vs High : p = %.4f\n', p23);

            posthoc_surface = table( ...
                {'Plat vs Medium'; 'Plat vs High'; 'Medium vs High'}, ...
                [p12; p13; p23], ...
                'VariableNames', {'Comparison', 'pValue'});
        end

        % Kruskal-Wallis (effet Groupe - entre-sujets)
        meanAcross = mean(allData, 2);
        [p_kruskal, ~, stats_kw] = kruskalwallis(meanAcross, groupLabels, 'off');
        fprintf('Kruskal-Wallis (Groupe) : p = %.4f\n', p_kruskal);

        posthoc_groupe = [];
        if p_kruskal < 0.05
            c_groupe_np = multcompare(stats_kw, 'Display', 'off');
            fprintf('Post-hoc Kruskal-Wallis :\n');
            disp(c_groupe_np);
            posthoc_groupe = c_groupe_np;
        end

        % Stockage
        resNP = struct();
        resNP.analysisType = 'nonparametric';
        resNP.dataWide     = T;
        resNP.meta.groupsOrder   = groupList;
        resNP.meta.surfacesOrder = Condition;
        resNP.descriptives = descriptives;

        resNP.surface.test    = 'Friedman';
        resNP.surface.p       = p_friedman;
        resNP.surface.posthoc = posthoc_surface;

        resNP.group.test      = 'KruskalWallis';
        resNP.group.p         = p_kruskal;
        resNP.group.posthoc   = posthoc_groupe;

        anovaResults.(matlab.lang.makeValidName(varName)) = resNP;
        continue;
    end

    % =========================================================
    % ===============  CAS PARAMÉTRIQUE  ======================
    % =========================================================
    try
        fprintf('➡️ Analyse PARAMÉTRIQUE (ANOVA mixte)\n');
        
        % Design intra-sujet
        within = table(categorical({'Plat'; 'Medium'; 'High'}), 'VariableNames', {'Surface'});

        % Modèle à mesures répétées
        rm = fitrm(T, 'Plat-High ~ Groupe', 'WithinDesign', within);

        % ANOVA intra-sujet (Surface + Interaction)
        ranovatbl = ranova(rm);
        fprintf('\n--- Table ANOVA à mesures répétées ---\n');
        disp(ranovatbl);
        
        % Test de sphéricité
        statsSphericity = mauchly(rm);
        fprintf('\n--- Test de sphéricité de Mauchly ---\n');
        disp(statsSphericity);

        % CORRECTION : ANOVA entre-sujets (Groupe)
        betweenTbl = anova(rm);
        fprintf('\n--- Table ANOVA entre-sujets ---\n');
        disp(betweenTbl);

        % CALCUL DES η² PARTIELS
        
        % --- η² pour Surface (intra-sujet) ---
idxSurf = strcmp(ranovatbl.Properties.RowNames, '(Intercept):Surface');
if any(idxSurf)
    SS_surf = ranovatbl.SumSq(idxSurf);
    idxErrSurf = find(idxSurf) + 1;
    if idxErrSurf <= height(ranovatbl)
        SS_err_surf = ranovatbl.SumSq(idxErrSurf);
        eta2_surface = SS_surf / (SS_surf + SS_err_surf);
    else
        eta2_surface = NaN;
    end
    pSurface = ranovatbl.pValueGG(idxSurf);
else
    eta2_surface = NaN;
    pSurface = NaN;
end

% --- η² pour Interaction ---
idxInt = strcmp(ranovatbl.Properties.RowNames, 'Groupe:Surface');
if any(idxInt)
    SS_int = ranovatbl.SumSq(idxInt);
    idxErrInt = find(idxInt) + 1;
    if idxErrInt <= height(ranovatbl)
        SS_err_int = ranovatbl.SumSq(idxErrInt);
        eta2_interaction = SS_int / (SS_int + SS_err_int);
    else
        eta2_interaction = NaN;
    end
    pInteraction = ranovatbl.pValueGG(idxInt);
else
    eta2_interaction = NaN;
    pInteraction = NaN;
end

% --- η² pour Groupe (entre-sujets) ---
pGroup = NaN;
eta2_group = NaN;

% La colonne s'appelle 'Between'
if ismember('Between', betweenTbl.Properties.VariableNames)
    % CORRECTION : Between est categorical, pas cell
    % Conversion en string pour comparaison
    betweenValues = string(betweenTbl.Between);
    
    % Recherche des lignes (case-insensitive)
    idxG = strcmpi(betweenValues, 'Groupe');
    idxErr = strcmpi(betweenValues, 'Error');
    
    if any(idxG) && any(idxErr)
        SS_group = betweenTbl.SumSq(idxG);
        SS_err_group = betweenTbl.SumSq(idxErr);
        eta2_group = SS_group / (SS_group + SS_err_group);
        pGroup = betweenTbl.pValue(idxG);
    else
        warning('Lignes Groupe ou Error introuvables dans betweenTbl pour %s', varName);
    end
elseif ismember('Source', betweenTbl.Properties.VariableNames)
    % Anciennes versions MATLAB
    sourceValues = string(betweenTbl.Source);
    idxG = strcmpi(sourceValues, 'Groupe');
    idxErr = strcmpi(sourceValues, 'Error');
    
    if any(idxG) && any(idxErr)
        SS_group = betweenTbl.SumSq(idxG);
        SS_err_group = betweenTbl.SumSq(idxErr);
        eta2_group = SS_group / (SS_group + SS_err_group);
        pGroup = betweenTbl.pValue(idxG);
    end
else
    % Fallback : extraction par indices de ligne
    if height(betweenTbl) >= 3
        SS_group = betweenTbl.SumSq(2);      % Ligne 2 = Groupe
        SS_err_group = betweenTbl.SumSq(3);  % Ligne 3 = Error
        pGroup = betweenTbl.pValue(2);
        eta2_group = SS_group / (SS_group + SS_err_group);
    else
        warning('Structure betweenTbl inconnue pour %s', varName);
    end
end

        fprintf('\n--- Tailles d''effet (η² partiel) ---\n');
        fprintf('Groupe     : η² = %.4f (p = %.4f)\n', eta2_group, pGroup);
        fprintf('Surface    : η² = %.4f (p = %.4f)\n', eta2_surface, pSurface);
        fprintf('Interaction: η² = %.4f (p = %.4f)\n', eta2_interaction, pInteraction);

        % POST-HOCS
        
        % Post-hoc Surface
        posthoc_surface = [];
        if pSurface < 0.05
            fprintf('\n--- Post-hoc Surface (comparaisons par paires) ---\n');
            posthoc_surface = multcompare(rm, 'Surface');
            disp(posthoc_surface);
        end

        % Post-hoc Interaction
        posthoc_interaction = [];
        if pInteraction < 0.05
            fprintf('\n--- Post-hoc Interaction (Surface × Groupe) ---\n');
            posthoc_interaction = multcompare(rm, 'Surface', 'By', 'Groupe');
            disp(posthoc_interaction);
        end

        % Post-hoc Groupe
        posthoc_groupe = [];
        if pGroup < 0.05
            fprintf('\n--- Post-hoc Groupe (comparaisons par paires) ---\n');
            meanAcross = mean(T{:, {'Plat','Medium','High'}}, 2);
            [~, ~, statsGroup] = anova1(meanAcross, T.Groupe, 'off');
            posthoc_groupe = multcompare(statsGroup, 'Display', 'off');
            disp(posthoc_groupe);
        end

        % STOCKAGE STRUCTURÉ
        res = struct();
        res.analysisType = 'parametric';
        res.dataWide     = T;
        res.meta.groupsOrder   = groupList;
        res.meta.surfacesOrder = Condition;
        res.descriptives = descriptives;
        res.ranova       = ranovatbl;
        res.rm           = rm;
        res.sphericity   = statsSphericity;
        res.between      = betweenTbl;

        % Surface
        res.surface.test    = 'RM-ANOVA';
        res.surface.p       = pSurface;
        res.surface.eta2    = eta2_surface;
        res.surface.posthoc = posthoc_surface;

        % Interaction
        res.interaction.test    = 'RM-ANOVA';
        res.interaction.p       = pInteraction;
        res.interaction.eta2    = eta2_interaction;
        res.interaction.posthoc = posthoc_interaction;

        % Groupe
        res.group.test    = 'Between-subjects ANOVA';
        res.group.p       = pGroup;
        res.group.eta2    = eta2_group;
        res.group.posthoc = posthoc_groupe;

        anovaResults.(matlab.lang.makeValidName(varName)) = res;

    catch ME
        fprintf('⚠️ Erreur pour %s : %s\n', varName, ME.message);
        fprintf('   Stack: %s\n', ME.stack(1).name);
        continue;
    end

end % fin boucle variables

%% SAUVEGARDE
save(fullfile(save_path, 'ANOVA_Results.mat'), 'anovaResults');
fprintf('\n✅ Résultats sauvegardés : %s\n', fullfile(save_path, 'ANOVA_Results.mat'));

%% EXPORT CSV RÉCAPITULATIF
fprintf('\n=== GÉNÉRATION DU TABLEAU RÉCAPITULATIF ===\n');

varNames = fieldnames(anovaResults);
recapTable = table();

for i = 1:length(varNames)
    varName = varNames{i};
    
    % Initialisation
    testUsed = 'Unknown';
    pGroup = NaN; eta2Group = NaN; starsGroup = '';
    pSurface = NaN; eta2Surface = NaN; starsSurface = '';
    pInteraction = NaN; eta2Interaction = NaN; starsInteraction = '';

    resVar = anovaResults.(varName);

    if isfield(resVar, 'analysisType')
        if strcmp(resVar.analysisType, 'parametric')
            testUsed = 'Parametric (Mixed ANOVA)';
            
            % Groupe
            if isfield(resVar, 'group')
                pGroup = resVar.group.p;
                eta2Group = resVar.group.eta2;
                starsGroup = repmat('*', 1, sum(pGroup < [0.05 0.01 0.001]));
            end
            
            % Surface
            if isfield(resVar, 'surface')
                pSurface = resVar.surface.p;
                eta2Surface = resVar.surface.eta2;
                starsSurface = repmat('*', 1, sum(pSurface < [0.05 0.01 0.001]));
            end
            
            % Interaction
            if isfield(resVar, 'interaction')
                pInteraction = resVar.interaction.p;
                eta2Interaction = resVar.interaction.eta2;
                starsInteraction = repmat('*', 1, sum(pInteraction < [0.05 0.01 0.001]));
            end
            
        elseif strcmp(resVar.analysisType, 'nonparametric')
            testUsed = 'Non-parametric';
            
            % Groupe
            if isfield(resVar, 'group')
                pGroup = resVar.group.p;
                starsGroup = repmat('*', 1, sum(pGroup < [0.05 0.01 0.001]));
            end
            
            % Surface
            if isfield(resVar, 'surface')
                pSurface = resVar.surface.p;
                starsSurface = repmat('*', 1, sum(pSurface < [0.05 0.01 0.001]));
            end
        end
    end

    % Ajout au tableau récapitulatif
    recapTable = [recapTable; table( ...
        {strrep(varName,'_',' ')}, ...
        {testUsed}, ...
        pGroup, {starsGroup}, eta2Group, ...
        pSurface, {starsSurface}, eta2Surface, ...
        pInteraction, {starsInteraction}, eta2Interaction, ...
        'VariableNames', {'Variable', 'Test', ...
                          'p_Groupe', 'Stars_Groupe', 'Eta2_Groupe', ...
                          'p_Surface', 'Stars_Surface', 'Eta2_Surface', ...
                          'p_Interaction', 'Stars_Interaction', 'Eta2_Interaction'})]; %#ok<AGROW>
end

% Export CSV
outputDir = fullfile(save_path, 'Anova_Recap');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

csvFilePath = fullfile(outputDir, 'ANOVA_Recap_Global.csv');
writetable(recapTable, csvFilePath);

fprintf('\n=== RÉSUMÉ FINAL ===\n');
fprintf('✅ Analyses ANOVA réalisées pour %d variables\n', length(varNames));
fprintf('✅ Tableau récapitulatif exporté : %s\n', csvFilePath);
fprintf('\n📊 Variables significatives (p < 0.05) :\n');

% Affichage des résultats significatifs
for i = 1:height(recapTable)
    if recapTable.p_Groupe(i) < 0.05 || recapTable.p_Surface(i) < 0.05 || recapTable.p_Interaction(i) < 0.05
        fprintf('  • %s\n', recapTable.Variable{i});
        if recapTable.p_Groupe(i) < 0.05
            fprintf('    - Groupe: p=%.4f, η²=%.3f %s\n', ...
                recapTable.p_Groupe(i), recapTable.Eta2_Groupe(i), recapTable.Stars_Groupe{i});
        end
        if recapTable.p_Surface(i) < 0.05
            fprintf('    - Surface: p=%.4f, η²=%.3f %s\n', ...
                recapTable.p_Surface(i), recapTable.Eta2_Surface(i), recapTable.Stars_Surface{i});
        end
        if recapTable.p_Interaction(i) < 0.05
            fprintf('    - Interaction: p=%.4f, η²=%.3f %s\n', ...
                recapTable.p_Interaction(i), recapTable.Eta2_Interaction(i), recapTable.Stars_Interaction{i});
        end
    end
end

disp('✅ Analyse terminée !');