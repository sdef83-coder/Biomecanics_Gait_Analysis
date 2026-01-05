%% LMM SUR LES PARAMÈTRES SPATIO-TEMPOREL - VERSION OPTIMISÉE
% Facteur 1: Groupe d'âge (between-subject)
% Facteur 2: Surface (within-subject, via (1|Participant))

clc; clear; close all;

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));

save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Statistics_LMM';
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

% Groupes et surfaces
groupList = {'JeunesEnfants', 'Enfants', 'Adolescents', 'Adultes'};
Condition = {'Plat', 'Medium', 'High'};

% Structures de sortie
anovaResultsLMM = struct();
recapTable = table();
posthocLMM = table();

alpha = 0.05;

%% BOUCLE PRINCIPALE POUR CHAQUE VARIABLE
for iVar = 1:length(variablesToAnalyze)
    varName = variablesToAnalyze{iVar};
    fprintf('\n=== LMM POUR : %s ===\n', varName);

    % ==============================
    % 1. Construire la table wide
    % ==============================
    allData = [];
    participantIDs = {};
    groupLabels = {};
    participantCounter = 1;

    for g = 1:length(groupList)
        groupName = groupList{g};

        if isfield(SpatioTemporalDATA, groupName) && ...
           isfield(SpatioTemporalDATA.(groupName), 'Plat')

            groupTable = SpatioTemporalDATA.(groupName).Plat;
            if ~isempty(groupTable) && any(strcmp(groupTable.Properties.VariableNames, varName))
                participants = unique(groupTable.Participant);

                for p = 1:length(participants)
                    participantID = participants{p};
                    participantData = [];
                    validData = true;

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
                                        validData = false; break;
                                    end
                                else
                                    validData = false; break;
                                end
                            else
                                validData = false; break;
                            end
                        else
                            validData = false; break;
                        end
                    end

                    if validData && length(participantData) == 3
                        allData = [allData; participantData]; %#ok<AGROW>
                        participantIDs{participantCounter,1} = participantID;
                        groupLabels{participantCounter,1} = groupName;
                        participantCounter = participantCounter + 1;
                    end
                end
            end
        end
    end

    nSubj = size(allData, 1);
    fprintf('Participants avec données complètes : %d\n', nSubj);

    if nSubj == 0
        warning('Aucune donnée complète pour %s → variable ignorée.', varName);
        continue;
    end

    groupLabelsCat = categorical(groupLabels(:), groupList, 'Ordinal', true);
    T_wide = table( ...
        categorical(participantIDs(:)), ...
        groupLabelsCat, ...
        allData(:,1), ...
        allData(:,2), ...
        allData(:,3), ...
        'VariableNames', {'Participant', 'Groupe', 'Plat', 'Medium', 'High'});

    % ==============================
    % 2. Construire le format long
    % ==============================
    T_long = table();
    for iC = 1:length(Condition)
        condName = Condition{iC};

        tmp = table();
        tmp.Participant = T_wide.Participant;
        tmp.Groupe      = T_wide.Groupe;
        tmp.Surface     = categorical( ...
                               repmat({condName}, nSubj, 1), ...
                               Condition, 'Ordinal', true);
        tmp.Y           = T_wide.(condName);

        T_long = [T_long; tmp]; %#ok<AGROW>
    end

    fprintf('Nombre total d''observations : %d\n', height(T_long));

    % ==============================
    % 3. Ajuster le LMM
    % ==============================
    try
        lme = fitlme(T_long, 'Y ~ Groupe*Surface + (1|Participant)');
        aovTbl = anova(lme, 'DFMethod', 'Satterthwaite');

        disp(aovTbl);

        % Extraire p-values et F
        pGroup = NaN; pSurface = NaN; pInteraction = NaN;
        FGroup = NaN; FSurface = NaN; FInteraction = NaN;

        if any(strcmp(aovTbl.Term, 'Groupe'))
            idx = strcmp(aovTbl.Term, 'Groupe');
            pGroup = aovTbl.pValue(idx);
            FGroup = aovTbl.FStat(idx);
        end
        if any(strcmp(aovTbl.Term, 'Surface'))
            idx = strcmp(aovTbl.Term, 'Surface');
            pSurface = aovTbl.pValue(idx);
            FSurface = aovTbl.FStat(idx);
        end
        if any(strcmp(aovTbl.Term, 'Groupe:Surface'))
            idx = strcmp(aovTbl.Term, 'Groupe:Surface');
            pInteraction = aovTbl.pValue(idx);
            FInteraction = aovTbl.FStat(idx);
        end

        fprintf('\n--- Résumé des effets fixes ---\n');
        fprintf('Groupe       : p = %.4f (F = %.3f)\n', pGroup, FGroup);
        fprintf('Surface      : p = %.4f (F = %.3f)\n', pSurface, FSurface);
        fprintf('Interaction  : p = %.4f (F = %.3f)\n', pInteraction, FInteraction);

        % Stocker
        res = struct();
        res.lme   = lme;
        res.aov   = aovTbl;
        res.pGroup = pGroup;
        res.pSurface = pSurface;
        res.pInteraction = pInteraction;
        res.FGroup = FGroup;
        res.FSurface = FSurface;
        res.FInteraction = FInteraction;

        anovaResultsLMM.(matlab.lang.makeValidName(varName)) = res;

        % Ajout au tableau récap
        recapTable = [recapTable; table( ...
            {strrep(varName,'_',' ')}, ...
            pGroup, FGroup, ...
            pSurface, FSurface, ...
            pInteraction, FInteraction, ...
            'VariableNames', {'Variable','p_Groupe','F_Groupe', ...
                              'p_Surface','F_Surface', ...
                              'p_Interaction','F_Interaction'})]; %#ok<AGROW>

        % ==============================
        % 4. POST-HOCS OPTIMISÉS
        % ==============================
        
        groupLevels = categories(T_long.Groupe);
        surfaceLevels = categories(T_long.Surface);
        nGroups = length(groupLevels);
        nSurfaces = length(surfaceLevels);
        
        coefNames = lme.CoefficientNames;
        df_residual = lme.DFE;
        
        fprintf('\n📊 Modèle : %d coefficients, DF = %.1f\n', ...
            lme.NumCoefficients, df_residual);

        % POST-HOC GROUPE
        if ~isnan(pGroup) && pGroup < alpha
            fprintf('\n   → Post-hoc Groupe\n');
            try
                nComp = (nGroups * (nGroups - 1)) / 2;
                pvals_raw = zeros(nComp, 1);
                estimates = zeros(nComp, 1);
                lowerCIs = zeros(nComp, 1);
                upperCIs = zeros(nComp, 1);
                compLabels = cell(nComp, 2);
                compIdx = 1;
                
                for i = 1:nGroups-1
                    for j = i+1:nGroups
                        g1 = groupLevels{i};
                        g2 = groupLevels{j};
                        
                        L = zeros(1, lme.NumCoefficients);
                        
                        if i == 1
                            idx2 = find(contains(coefNames, ['Groupe_' char(g2)]), 1);
                            if ~isempty(idx2)
                                L(idx2) = 1;
                            else
                                compIdx = compIdx + 1;
                                continue;
                            end
                        else
                            idx1 = find(contains(coefNames, ['Groupe_' char(g1)]), 1);
                            idx2 = find(contains(coefNames, ['Groupe_' char(g2)]), 1);
                            
                            if ~isempty(idx1) && ~isempty(idx2)
                                L(idx1) = 1;
                                L(idx2) = -1;
                            else
                                compIdx = compIdx + 1;
                                continue;
                            end
                        end
                        
                        [p, est, lCI, uCI, ~] = testContrast(lme, L);
                        
                        pvals_raw(compIdx) = p;
                        estimates(compIdx) = est;
                        lowerCIs(compIdx) = lCI;
                        upperCIs(compIdx) = uCI;
                        compLabels{compIdx, 1} = string(g1);
                        compLabels{compIdx, 2} = string(g2);
                        
                        compIdx = compIdx + 1;
                    end
                end
                
                % Supprimer entrées vides
                validIdx = pvals_raw > 0;
                pvals_raw = pvals_raw(validIdx);
                estimates = estimates(validIdx);
                lowerCIs = lowerCIs(validIdx);
                upperCIs = upperCIs(validIdx);
                compLabels = compLabels(validIdx, :);
                
                % Correction de Holm
                [sortedP, idxSort] = sort(pvals_raw);
                m = length(pvals_raw);
                pvals_adj = zeros(m, 1);
                
                for k = 1:m
                    pvals_adj(idxSort(k)) = min(sortedP(k) * (m - k + 1), 1);
                end
                
                for k = m-1:-1:1
                    pvals_adj(idxSort(k)) = min(pvals_adj(idxSort(k)), pvals_adj(idxSort(k+1)));
                end
                
                % Stocker
                for k = 1:m
                    newRow = table( ...
                        string(varName), ...
                        "Groupe", ...
                        compLabels{k,1}, compLabels{k,2}, ...
                        estimates(k), ...
                        lowerCIs(k), upperCIs(k), ...
                        pvals_raw(k), ...
                        pvals_adj(k), ...
                        "LMM_Holm", ...
                        'VariableNames', {'Variable','Effect','Level1','Level2', ...
                                          'Estimate','LowerCI','UpperCI','pValue_raw','pValue_adj','Method'});
                    posthocLMM = [posthocLMM; newRow]; %#ok<AGROW>
                    
                    if pvals_adj(k) < 0.05
                        fprintf('    %s vs %s: p_adj = %.4f *\n', ...
                            compLabels{k,1}, compLabels{k,2}, pvals_adj(k));
                    end
                end
                
            catch ME
                fprintf('⚠️ Erreur post-hoc Groupe : %s\n', ME.message);
            end
        end

        % POST-HOC SURFACE
        if ~isnan(pSurface) && pSurface < alpha
            fprintf('\n   → Post-hoc Surface\n');
            try
                nComp = (nSurfaces * (nSurfaces - 1)) / 2;
                pvals_raw = zeros(nComp, 1);
                estimates = zeros(nComp, 1);
                lowerCIs = zeros(nComp, 1);
                upperCIs = zeros(nComp, 1);
                compLabels = cell(nComp, 2);
                compIdx = 1;
                
                for i = 1:nSurfaces-1
                    for j = i+1:nSurfaces
                        s1 = surfaceLevels{i};
                        s2 = surfaceLevels{j};
                        
                        L = zeros(1, lme.NumCoefficients);
                        
                        if i == 1
                            idx2 = find(contains(coefNames, ['Surface_' char(s2)]), 1);
                            if ~isempty(idx2)
                                L(idx2) = 1;
                            else
                                compIdx = compIdx + 1;
                                continue;
                            end
                        else
                            idx1 = find(contains(coefNames, ['Surface_' char(s1)]), 1);
                            idx2 = find(contains(coefNames, ['Surface_' char(s2)]), 1);
                            
                            if ~isempty(idx1) && ~isempty(idx2)
                                L(idx1) = 1;
                                L(idx2) = -1;
                            else
                                compIdx = compIdx + 1;
                                continue;
                            end
                        end
                        
                        [p, est, lCI, uCI, ~] = testContrast(lme, L);
                        
                        pvals_raw(compIdx) = p;
                        estimates(compIdx) = est;
                        lowerCIs(compIdx) = lCI;
                        upperCIs(compIdx) = uCI;
                        compLabels{compIdx, 1} = string(s1);
                        compLabels{compIdx, 2} = string(s2);
                        
                        compIdx = compIdx + 1;
                    end
                end
                
                validIdx = pvals_raw > 0;
                pvals_raw = pvals_raw(validIdx);
                estimates = estimates(validIdx);
                lowerCIs = lowerCIs(validIdx);
                upperCIs = upperCIs(validIdx);
                compLabels = compLabels(validIdx, :);
                
                [sortedP, idxSort] = sort(pvals_raw);
                m = length(pvals_raw);
                pvals_adj = zeros(m, 1);
                
                for k = 1:m
                    pvals_adj(idxSort(k)) = min(sortedP(k) * (m - k + 1), 1);
                end
                
                for k = m-1:-1:1
                    pvals_adj(idxSort(k)) = min(pvals_adj(idxSort(k)), pvals_adj(idxSort(k+1)));
                end
                
                for k = 1:m
                    newRow = table( ...
                        string(varName), ...
                        "Surface", ...
                        compLabels{k,1}, compLabels{k,2}, ...
                        estimates(k), ...
                        lowerCIs(k), upperCIs(k), ...
                        pvals_raw(k), ...
                        pvals_adj(k), ...
                        "LMM_Holm", ...
                        'VariableNames', {'Variable','Effect','Level1','Level2', ...
                                          'Estimate','LowerCI','UpperCI','pValue_raw','pValue_adj','Method'});
                    posthocLMM = [posthocLMM; newRow]; %#ok<AGROW>
                    
                    if pvals_adj(k) < 0.05
                        fprintf('    %s vs %s: p_adj = %.4f *\n', ...
                            compLabels{k,1}, compLabels{k,2}, pvals_adj(k));
                    end
                end
                
            catch ME
                fprintf('⚠️ Erreur post-hoc Surface : %s\n', ME.message);
            end
        end

        % POST-HOC INTERACTION
        if ~isnan(pInteraction) && pInteraction < alpha
            fprintf('\n   → Post-hoc Interaction\n');
            try
                for g = 1:nGroups
                    gName = groupLevels{g};
                    fprintf('    Groupe: %s\n', gName);
                    
                    nComp_int = (nSurfaces * (nSurfaces - 1)) / 2;
                    pvals_raw_g = zeros(nComp_int, 1);
                    estimates_g = zeros(nComp_int, 1);
                    lowerCIs_g = zeros(nComp_int, 1);
                    upperCIs_g = zeros(nComp_int, 1);
                    compLabels_g = cell(nComp_int, 2);
                    compIdx_g = 1;
                    
                    for i = 1:nSurfaces-1
                        for j = i+1:nSurfaces
                            s1 = surfaceLevels{i};
                            s2 = surfaceLevels{j};
                            
                            L = zeros(1, lme.NumCoefficients);
                            
                            if g == 1 && i == 1
                                idx_surf = find(contains(coefNames, ['Surface_' char(s2)]), 1);
                                if ~isempty(idx_surf)
                                    L(idx_surf) = 1;
                                else
                                    compIdx_g = compIdx_g + 1;
                                    continue;
                                end
                            elseif g == 1
                                idx_s1 = find(contains(coefNames, ['Surface_' char(s1)]), 1);
                                idx_s2 = find(contains(coefNames, ['Surface_' char(s2)]), 1);
                                if ~isempty(idx_s1) && ~isempty(idx_s2)
                                    L(idx_s1) = 1;
                                    L(idx_s2) = -1;
                                else
                                    compIdx_g = compIdx_g + 1;
                                    continue;
                                end
                            else
                                if i == 1
                                    idx_surf = find(contains(coefNames, ['Surface_' char(s2)]), 1);
                                    idx_int = find(contains(coefNames, ['Groupe_' char(gName)]) & ...
                                                  contains(coefNames, ['Surface_' char(s2)]), 1);
                                    
                                    if ~isempty(idx_surf), L(idx_surf) = 1; end
                                    if ~isempty(idx_int), L(idx_int) = 1; end
                                    
                                    if isempty(idx_surf) && isempty(idx_int)
                                        compIdx_g = compIdx_g + 1;
                                        continue;
                                    end
                                else
                                    idx_s1 = find(contains(coefNames, ['Surface_' char(s1)]), 1);
                                    idx_s2 = find(contains(coefNames, ['Surface_' char(s2)]), 1);
                                    idx_int_s1 = find(contains(coefNames, ['Groupe_' char(gName)]) & ...
                                                     contains(coefNames, ['Surface_' char(s1)]), 1);
                                    idx_int_s2 = find(contains(coefNames, ['Groupe_' char(gName)]) & ...
                                                     contains(coefNames, ['Surface_' char(s2)]), 1);
                                    
                                    if ~isempty(idx_s1), L(idx_s1) = 1; end
                                    if ~isempty(idx_s2), L(idx_s2) = -1; end
                                    if ~isempty(idx_int_s1), L(idx_int_s1) = 1; end
                                    if ~isempty(idx_int_s2), L(idx_int_s2) = -1; end
                                    
                                    if all(L == 0)
                                        compIdx_g = compIdx_g + 1;
                                        continue;
                                    end
                                end
                            end
                            
                            [p, est, lCI, uCI, ~] = testContrast(lme, L);
                            
                            pvals_raw_g(compIdx_g) = p;
                            estimates_g(compIdx_g) = est;
                            lowerCIs_g(compIdx_g) = lCI;
                            upperCIs_g(compIdx_g) = uCI;
                            compLabels_g{compIdx_g, 1} = string(s1);
                            compLabels_g{compIdx_g, 2} = string(s2);
                            
                            compIdx_g = compIdx_g + 1;
                        end
                    end
                    
                    validIdx_g = pvals_raw_g > 0;
                    pvals_raw_g = pvals_raw_g(validIdx_g);
                    estimates_g = estimates_g(validIdx_g);
                    lowerCIs_g = lowerCIs_g(validIdx_g);
                    upperCIs_g = upperCIs_g(validIdx_g);
                    compLabels_g = compLabels_g(validIdx_g, :);
                    
                    if isempty(pvals_raw_g)
                        fprintf('      → Aucun contraste valide\n');
                        continue;
                    end
                    
                    [sortedP_g, idxSort_g] = sort(pvals_raw_g);
                    m_g = length(pvals_raw_g);
                    pvals_adj_g = zeros(m_g, 1);
                    
                    for k = 1:m_g
                        pvals_adj_g(idxSort_g(k)) = min(sortedP_g(k) * (m_g - k + 1), 1);
                    end
                    
                    for k = m_g-1:-1:1
                        pvals_adj_g(idxSort_g(k)) = min(pvals_adj_g(idxSort_g(k)), ...
                                                         pvals_adj_g(idxSort_g(k+1)));
                    end
                    
                    for k = 1:m_g
                        newRow = table( ...
                            string(varName), ...
                            "Surface|Groupe", ...
                            compLabels_g{k,1} + " | " + string(gName), ...
                            compLabels_g{k,2} + " | " + string(gName), ...
                            estimates_g(k), ...
                            lowerCIs_g(k), upperCIs_g(k), ...
                            pvals_raw_g(k), ...
                            pvals_adj_g(k), ...
                            "LMM_Contraste_Holm", ...
                            'VariableNames', {'Variable','Effect','Level1','Level2', ...
                                              'Estimate','LowerCI','UpperCI','pValue_raw','pValue_adj','Method'});
                        posthocLMM = [posthocLMM; newRow]; %#ok<AGROW>
                        
                        if pvals_adj_g(k) < 0.05
                            fprintf('      %s vs %s: p_adj = %.4f *\n', ...
                                compLabels_g{k,1}, compLabels_g{k,2}, pvals_adj_g(k));
                        end
                    end
                end
                
            catch ME
                fprintf('⚠️ Erreur post-hoc interaction : %s\n', ME.message);
            end
        end

    catch ME
        fprintf('⚠️ Erreur LMM pour %s : %s\n', varName, ME.message);
        if ~isempty(ME.stack)
            fprintf('   Stack: %s (ligne %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        continue;
    end
end

%% SAUVEGARDE DES RÉSULTATS LMM
save(fullfile(save_path, 'LMM_Results.mat'), 'anovaResultsLMM');
fprintf('\n✅ Résultats LMM sauvegardés : %s\n', fullfile(save_path, 'LMM_Results.mat'));

%% EXPORT CSV RÉCAPITULATIF DES EFFETS FIXES
outputDir = fullfile(save_path, 'LMM_Recap');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

csvRecapPath = fullfile(outputDir, 'LMM_Recap_FixedEffects.csv');

% indicateur de significativité dans le récapitulatif
stars = @(p) repmat('*',1,sum(p < [0.05 0.01 0.001]));

recapTable.Sig_Groupe      = cellfun(stars, num2cell(recapTable.p_Groupe), 'UniformOutput', false);
recapTable.Sig_Surface     = cellfun(stars, num2cell(recapTable.p_Surface), 'UniformOutput', false);
recapTable.Sig_Interaction = cellfun(stars, num2cell(recapTable.p_Interaction), 'UniformOutput', false);

writetable(recapTable, csvRecapPath);

fprintf('\n=== RÉSUMÉ FINAL LMM ===\n');
fprintf('✅ LMM réalisés pour %d variables\n', height(recapTable));
fprintf('✅ Tableau récapitulatif exporté : %s\n', csvRecapPath);

%% EXPORT CSV DES POST-HOCS
if ~isempty(posthocLMM)
    Stars = cell(height(posthocLMM),1);
    for i = 1:height(posthocLMM)
        % utiliser les p-values ajustées
        p = posthocLMM.pValue_adj(i);
        nStar = sum(p < [0.05 0.01 0.001]);
        Stars{i} = repmat('*',1,nStar);
    end
    posthocLMM.Stars = Stars;
end

csvPosthocPath = fullfile(outputDir, 'LMM_PostHoc_All.csv');
writetable(posthocLMM, csvPosthocPath);
fprintf('📄 Post-hocs LMM exportés dans : %s\n', csvPosthocPath);

%% Fonction helper
function [pval, est, lowerCI, upperCI, df] = testContrast(lme, L)
        % Test du contraste
        [pval, F, DF1, DF2] = coefTest(lme, L);
        
        % Estimation de la différence
        est = L * lme.Coefficients.Estimate;
        
        % Erreur standard
        C = lme.CoefficientCovariance;
        se = sqrt(L * C * L');
        
        % Distribution t (pas normale!)
        df = DF2;
        t_crit = tinv(0.975, df);
        
        % Intervalles de confiance
        lowerCI = est - t_crit * se;
        upperCI = est + t_crit * se;
    end