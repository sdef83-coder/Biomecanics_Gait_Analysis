%% LMM Surface x Groupe (standard littérature) - robuste aux tailles inégales et données manquantes
% Design:
%   - Groupe (between-subject) : JeunesEnfants, Enfants, Adolescents, Adultes
%   - Surface (within-subject) : Plat, Medium, High
% Modèle principal (standard) :
%   Y ~ Groupe*Surface + (1|Participant)
% Option recommandée :
%   comparer avec Y ~ Groupe*Surface + (1 + Surface|Participant)
% Post-hoc (standard) :
%   estimated marginal means (EMMs) + contrasts + correction 

clc; clear; close all;

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));

save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Statistics_LMM_STANDARD';
if ~exist(save_path, 'dir'); mkdir(save_path); end

load('SpatioTemporalDATA.mat');

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

groupList   = {'JeunesEnfants','Enfants','Adolescents','Adultes'};
surfaceList = {'Plat','Medium','High'};

alpha = 0.05;

resultsLMM = struct();
recapTable = table();
posthocAll = table();

for iVar = 1:numel(variablesToAnalyze)
    varName = variablesToAnalyze{iVar};
    fprintf('\n====================================================\n');
    fprintf('LMM (standard) pour : %s\n', varName);
    fprintf('====================================================\n');

    % =========================================================
    % 1) Construire la table LONG en gardant toutes les obs dispo
    % =========================================================
    T_long = buildLongTableFromStruct(SpatioTemporalDATA, varName, groupList, surfaceList);

    if isempty(T_long) || height(T_long) < 6
        warning('Pas assez de données pour %s → ignorée.', varName);
        continue;
    end

    % Garde-fous simples (standard)
    nSubj = numel(categories(T_long.Participant));
    if nSubj < 6
        warning('Trop peu de participants (%d) pour %s → ignorée.', nSubj, varName);
        continue;
    end

    fprintf('Participants uniques : %d\n', nSubj);
    fprintf('Observations totales : %d\n', height(T_long));

    % =========================================================
    % 2) Ajuster LMM - modèle standard + option random slope
    % =========================================================
    try
        % (A) Modèle standard : random intercept
        lme_RI = fitlme(T_long, 'Y ~ Groupe*Surface + (1|Participant)', ...
            'FitMethod','REML');

        % (B) Option recommandée : random slope de Surface
        % Pour comparer des structures de random effects, on compare en ML.
        lme_RI_ML = fitlme(T_long, 'Y ~ Groupe*Surface + (1|Participant)', ...
            'FitMethod','ML');
        lme_RS_ML = [];
        bestModel = lme_RI;         % par défaut
        bestModelName = "RI_REML";  % random intercept

        try
            lme_RS_ML = fitlme(T_long, 'Y ~ Groupe*Surface + (1 + Surface|Participant)', ...
                'FitMethod','ML');

            % Sélection simple (standard) : AIC plus faible en ML
            if lme_RS_ML.ModelCriterion.AIC + 2 < lme_RI_ML.ModelCriterion.AIC
                % Refit en REML pour estimation finale (standard)
                bestModel = fitlme(T_long, 'Y ~ Groupe*Surface + (1 + Surface|Participant)', ...
                    'FitMethod','REML');
                bestModelName = "RS_REML";
            end
        catch
            % si le modèle slope ne converge pas, on reste sur RI
        end

        lme = bestModel;

        fprintf('Modèle retenu : %s\n', bestModelName);
        fprintf('Formule : %s\n', lme.Formula);

        % =========================================================
        % 3) Table ANOVA des effets fixes (Satterthwaite)
        % =========================================================
        aov = anova(lme, 'DFMethod','Satterthwaite');
        disp(aov);

        [pGroup, FGroup]         = getTermStats_STD(aov, 'Groupe');
        [pSurface, FSurface]     = getTermStats_STD(aov, 'Surface');
        [pInter, FInter]         = getTermStats_STD(aov, 'Groupe:Surface');

        fprintf('\nEffets fixes (Satterthwaite):\n');
        fprintf('  Groupe       : p=%.4f, F=%.3f\n', pGroup,   FGroup);
        fprintf('  Surface      : p=%.4f, F=%.3f\n', pSurface, FSurface);
        fprintf('  Interaction  : p=%.4f, F=%.3f\n', pInter,   FInter);

% =========================================================
% POST-HOC (standard)
%   - If interaction is significant:
%         (A) Surface within each Groupe (Holm within-group)
%         (B) Groupe within each Surface (Holm within-surface)
%   - Else:
%         Post-hocs main effects (Groupe / Surface) ONLY if significant
% =========================================================

% Fixed-effects components (used everywhere below)
D      = designMatrix(lme,'Fixed');         % Nobs x P
beta   = fixedEffects(lme);                 % P x 1
C      = lme.CoefficientCovariance;         % P x P
df_res = lme.DFE;

groupLevels   = categories(T_long.Groupe);
surfaceLevels = categories(T_long.Surface);

if ~isnan(pInter) && pInter < alpha

    % =========================================================
    % (A) POST-HOC INTERACTION: Surface within each Groupe
    % =========================================================
    fprintf('\n   → Post-hoc Interaction A: Surface within each Groupe (Holm)\n');

    for g = 1:numel(groupLevels)
        gName = groupLevels{g};
        fprintf('    Groupe: %s\n', gName);

        % xbar: one "design row" per surface within this group
        xbar = NaN(numel(surfaceLevels), size(D,2));
        for s = 1:numel(surfaceLevels)
            sName = surfaceLevels{s};
            idx = (T_long.Groupe == gName) & (T_long.Surface == sName);
            if any(idx)
                xbar(s,:) = mean(D(idx,:), 1);
            end
        end

        % pairwise surfaces
        pairs   = nchoosek(1:numel(surfaceLevels), 2);
        compLab = strings(size(pairs,1),1);
        p_raw   = NaN(size(pairs,1),1);
        est     = NaN(size(pairs,1),1);
        lCI     = NaN(size(pairs,1),1);
        uCI     = NaN(size(pairs,1),1);

        for k = 1:size(pairs,1)
            i1 = pairs(k,1); i2 = pairs(k,2);
            s1 = string(surfaceLevels{i1}); s2 = string(surfaceLevels{i2});
            compLab(k) = s1 + " vs " + s2;

            if any(isnan(xbar(i1,:))) || any(isnan(xbar(i2,:)))
                continue;
            end

            L = xbar(i1,:) - xbar(i2,:);

            p_raw(k) = coefTest(lme, L);

            est(k) = L * beta;
            se     = sqrt(L * C * L');
            tcrit  = tinv(0.975, df_res);
            lCI(k) = est(k) - tcrit * se;
            uCI(k) = est(k) + tcrit * se;
        end

        valid = ~isnan(p_raw);
        if ~any(valid)
            fprintf('      → Aucun contraste valide\n');
            continue;
        end

        p_adj = NaN(size(p_raw));
        p_adj(valid) = holmAdjustVector(p_raw(valid));

        for k = 1:numel(p_raw)
            if isnan(p_raw(k)); continue; end

            newRow = table( ...
                string(varName), ...
                "Interaction_SurfaceWithinGroup", ...
                string(gName), ...
                compLab(k), ...
                est(k), lCI(k), uCI(k), ...
                p_raw(k), p_adj(k), ...
                "Holm_withinGroup", ...
                'VariableNames', {'Variable','Effect','Group','Comparison', ...
                                  'Estimate','LowerCI','UpperCI','pValue_raw','pValue_adj','Method'} );

            posthocAll = [posthocAll; newRow]; %#ok<AGROW>

            if p_adj(k) < alpha
                fprintf('      %s: p_adj=%.4f *  (Est=%.4g, CI=[%.4g; %.4g])\n', ...
                    compLab(k), p_adj(k), est(k), lCI(k), uCI(k));
            else
                fprintf('      %s: p_adj=%.4f    (Est=%.4g)\n', ...
                    compLab(k), p_adj(k), est(k));
            end
        end
    end

    % =========================================================
    % (B) POST-HOC INTERACTION: Groupe within each Surface
    % =========================================================
    fprintf('\n   → Post-hoc Interaction B: Groupe within each Surface (Holm)\n');

    for s = 1:numel(surfaceLevels)
        sName = surfaceLevels{s};
        fprintf('    Surface: %s\n', sName);

        % xbar: one "design row" per group within this surface
        xbar = NaN(numel(groupLevels), size(D,2));
        for g = 1:numel(groupLevels)
            gName = groupLevels{g};
            idx = (T_long.Groupe == gName) & (T_long.Surface == sName);
            if any(idx)
                xbar(g,:) = mean(D(idx,:), 1);
            end
        end

        % pairwise groups
        pairs   = nchoosek(1:numel(groupLevels), 2);
        compLab = strings(size(pairs,1),1);
        p_raw   = NaN(size(pairs,1),1);
        est     = NaN(size(pairs,1),1);
        lCI     = NaN(size(pairs,1),1);
        uCI     = NaN(size(pairs,1),1);

        for k = 1:size(pairs,1)
            i1 = pairs(k,1); i2 = pairs(k,2);
            g1 = string(groupLevels{i1}); g2 = string(groupLevels{i2});
            compLab(k) = g1 + " vs " + g2;

            if any(isnan(xbar(i1,:))) || any(isnan(xbar(i2,:)))
                continue;
            end

            L = xbar(i1,:) - xbar(i2,:);

            p_raw(k) = coefTest(lme, L);

            est(k) = L * beta;
            se     = sqrt(L * C * L');
            tcrit  = tinv(0.975, df_res);
            lCI(k) = est(k) - tcrit * se;
            uCI(k) = est(k) + tcrit * se;
        end

        valid = ~isnan(p_raw);
        if ~any(valid)
            fprintf('      → Aucun contraste valide\n');
            continue;
        end

        p_adj = NaN(size(p_raw));
        p_adj(valid) = holmAdjustVector(p_raw(valid));

        for k = 1:numel(p_raw)
            if isnan(p_raw(k)); continue; end

            newRow = table( ...
                string(varName), ...
                "Interaction_GroupWithinSurface", ...
                "All", ...                        % colonne Group: on met "All" car c'est l'effet "Groupe"
                "Surface=" + string(sName) + ": " + compLab(k), ...
                est(k), lCI(k), uCI(k), ...
                p_raw(k), p_adj(k), ...
                "Holm_withinSurface", ...
                'VariableNames', {'Variable','Effect','Group','Comparison', ...
                                  'Estimate','LowerCI','UpperCI','pValue_raw','pValue_adj','Method'} );

            posthocAll = [posthocAll; newRow]; %#ok<AGROW>

            if p_adj(k) < alpha
                fprintf('      %s: p_adj=%.4f *  (Est=%.4g, CI=[%.4g; %.4g])\n', ...
                    compLab(k), p_adj(k), est(k), lCI(k), uCI(k));
            else
                fprintf('      %s: p_adj=%.4f    (Est=%.4g)\n', ...
                    compLab(k), p_adj(k), est(k));
            end
        end
    end

else
    % =========================================================
    % NO interaction (or NaN): main effects post-hocs if significant
    % =========================================================

    % ---------- Post-hoc Groupe (moyenné sur Surface) ----------
    if ~isnan(pGroup) && pGroup < alpha
        fprintf('\n   → Post-hoc Groupe (moyenné sur Surface, Holm)\n');

        % design row marginale par groupe (moyenne égale sur surfaces, conservateur)
        xbarG = NaN(numel(groupLevels), size(D,2));

        for g = 1:numel(groupLevels)
            gName = groupLevels{g};

            rowsS = NaN(numel(surfaceLevels), size(D,2));
            for s = 1:numel(surfaceLevels)
                sName = surfaceLevels{s};
                idx = (T_long.Groupe == gName) & (T_long.Surface == sName);
                if any(idx)
                    rowsS(s,:) = mean(D(idx,:),1);
                end
            end

            % Conservateur: il faut que toutes les surfaces existent pour ce groupe
            if all(~any(isnan(rowsS),2))
                xbarG(g,:) = mean(rowsS, 1);
            end
        end

        pairs   = nchoosek(1:numel(groupLevels), 2);
        compLab = strings(size(pairs,1),1);
        p_raw   = NaN(size(pairs,1),1);
        est     = NaN(size(pairs,1),1);
        lCI     = NaN(size(pairs,1),1);
        uCI     = NaN(size(pairs,1),1);

        for k = 1:size(pairs,1)
            i1 = pairs(k,1); i2 = pairs(k,2);
            g1 = string(groupLevels{i1}); g2 = string(groupLevels{i2});
            compLab(k) = g1 + " vs " + g2;

            if any(isnan(xbarG(i1,:))) || any(isnan(xbarG(i2,:)))
                continue;
            end

            L = xbarG(i1,:) - xbarG(i2,:);

            p_raw(k) = coefTest(lme, L);

            est(k) = L * beta;
            se     = sqrt(L * C * L');
            tcrit  = tinv(0.975, df_res);
            lCI(k) = est(k) - tcrit * se;
            uCI(k) = est(k) + tcrit * se;
        end

        valid = ~isnan(p_raw);
        if any(valid)
            p_adj = NaN(size(p_raw));
            p_adj(valid) = holmAdjustVector(p_raw(valid));
        else
            p_adj = NaN(size(p_raw));
        end

        for k = 1:numel(p_raw)
            if isnan(p_raw(k)); continue; end

            newRow = table( ...
                string(varName), ...
                "MainEffect_Groupe", ...
                "All", ...
                compLab(k), ...
                est(k), lCI(k), uCI(k), ...
                p_raw(k), p_adj(k), ...
                "Holm", ...
                'VariableNames', {'Variable','Effect','Group','Comparison', ...
                                  'Estimate','LowerCI','UpperCI','pValue_raw','pValue_adj','Method'} );

            posthocAll = [posthocAll; newRow]; %#ok<AGROW>

            if p_adj(k) < alpha
                fprintf('      %s: p_adj=%.4f *  (Est=%.4g, CI=[%.4g; %.4g])\n', ...
                    compLab(k), p_adj(k), est(k), lCI(k), uCI(k));
            else
                fprintf('      %s: p_adj=%.4f    (Est=%.4g)\n', ...
                    compLab(k), p_adj(k), est(k));
            end
        end
    end

    % ---------- Post-hoc Surface (moyenné sur Groupe) ----------
    if ~isnan(pSurface) && pSurface < alpha
        fprintf('\n   → Post-hoc Surface (moyenné sur Groupe, Holm)\n');

        % design row marginale par surface (moyenne égale sur groupes, conservateur)
        xbarS = NaN(numel(surfaceLevels), size(D,2));

        for s = 1:numel(surfaceLevels)
            sName = surfaceLevels{s};

            rowsG = NaN(numel(groupLevels), size(D,2));
            for g = 1:numel(groupLevels)
                gName = groupLevels{g};
                idx = (T_long.Groupe == gName) & (T_long.Surface == sName);
                if any(idx)
                    rowsG(g,:) = mean(D(idx,:),1);
                end
            end

            % Conservateur: il faut que tous les groupes existent pour cette surface
            if all(~any(isnan(rowsG),2))
                xbarS(s,:) = mean(rowsG, 1);
            end
        end

        pairs   = nchoosek(1:numel(surfaceLevels), 2);
        compLab = strings(size(pairs,1),1);
        p_raw   = NaN(size(pairs,1),1);
        est     = NaN(size(pairs,1),1);
        lCI     = NaN(size(pairs,1),1);
        uCI     = NaN(size(pairs,1),1);

        for k = 1:size(pairs,1)
            i1 = pairs(k,1); i2 = pairs(k,2);
            s1 = string(surfaceLevels{i1}); s2 = string(surfaceLevels{i2});
            compLab(k) = s1 + " vs " + s2;

            if any(isnan(xbarS(i1,:))) || any(isnan(xbarS(i2,:)))
                continue;
            end

            L = xbarS(i1,:) - xbarS(i2,:);

            p_raw(k) = coefTest(lme, L);

            est(k) = L * beta;
            se     = sqrt(L * C * L');
            tcrit  = tinv(0.975, df_res);
            lCI(k) = est(k) - tcrit * se;
            uCI(k) = est(k) + tcrit * se;
        end

        valid = ~isnan(p_raw);
        if any(valid)
            p_adj = NaN(size(p_raw));
            p_adj(valid) = holmAdjustVector(p_raw(valid));
        else
            p_adj = NaN(size(p_raw));
        end

        for k = 1:numel(p_raw)
            if isnan(p_raw(k)); continue; end

            newRow = table( ...
                string(varName), ...
                "MainEffect_Surface", ...
                "All", ...
                compLab(k), ...
                est(k), lCI(k), uCI(k), ...
                p_raw(k), p_adj(k), ...
                "Holm", ...
                'VariableNames', {'Variable','Effect','Group','Comparison', ...
                                  'Estimate','LowerCI','UpperCI','pValue_raw','pValue_adj','Method'} );

            posthocAll = [posthocAll; newRow]; %#ok<AGROW>

            if p_adj(k) < alpha
                fprintf('      %s: p_adj=%.4f *  (Est=%.4g, CI=[%.4g; %.4g])\n', ...
                    compLab(k), p_adj(k), est(k), lCI(k), uCI(k));
            else
                fprintf('      %s: p_adj=%.4f    (Est=%.4g)\n', ...
                    compLab(k), p_adj(k), est(k));
            end
        end
    end
end

        % =========================================================
        % 5) Stockage + récapitulatif
        % =========================================================
        res = struct();
        res.varName = varName;
        res.modelName = bestModelName;
        res.lme = lme;
        res.aov = aov;
        res.pGroup = pGroup; res.FGroup = FGroup;
        res.pSurface = pSurface; res.FSurface = FSurface;
        res.pInteraction = pInter; res.FInteraction = FInter;
        %res.EMMs = grid;

        resultsLMM.(matlab.lang.makeValidName(varName)) = res;

        recapTable = [recapTable; table( ...
            string(varName), string(bestModelName), ...
            pGroup, FGroup, ...
            pSurface, FSurface, ...
            pInter, FInter, ...
            'VariableNames', {'Variable','Model', ...
                              'p_Groupe','F_Groupe', ...
                              'p_Surface','F_Surface', ...
                              'p_Interaction','F_Interaction'})]; %#ok<AGROW>

    catch ME
        fprintf('⚠️ Erreur LMM pour %s : %s\n', varName, ME.message);
        if ~isempty(ME.stack)
            fprintf('   -> %s (ligne %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        continue;
    end
end

%% Sauvegarde
save(fullfile(save_path, 'LMM_Results.mat'), 'resultsLMM');

%% Exports CSV
outDir = fullfile(save_path, 'Exports');
if ~exist(outDir,'dir'); mkdir(outDir); end

% Récap
recapTable.Sig_Groupe      = starsCol(recapTable.p_Groupe);
recapTable.Sig_Surface     = starsCol(recapTable.p_Surface);
recapTable.Sig_Interaction = starsCol(recapTable.p_Interaction);
writetable(recapTable, fullfile(outDir, 'LMM_Recap_FixedEffects.csv'));

% Post-hoc
if ~isempty(posthocAll)
    writetable(posthocAll, fullfile(outDir, 'LMM_PostHoc_EMM_Holm.csv'));
end

fprintf('\n✅ Terminé.\n');
fprintf('  - MAT: %s\n', fullfile(save_path, 'LMM_Results_STANDARD.mat'));
fprintf('  - CSV recap: %s\n', fullfile(outDir, 'LMM_Recap_FixedEffects.csv'));
if ~isempty(posthocAll)
    fprintf('  - CSV posthoc: %s\n', fullfile(outDir, 'LMM_PostHoc_EMM_Holm.csv'));
end

%% Fonctions locales (helpers)

function T_long = buildLongTableFromStruct(S, varName, groupList, surfaceList)
% Construit un format long en conservant toutes les observations disponibles.
% Hypothèse: S.(group).(surface) est une table contenant 'Participant' + varName.

    rows = {};
    r = 0;

    for g = 1:numel(groupList)
        groupName = groupList{g};
        if ~isfield(S, groupName); continue; end

        % On prend l'ensemble des participants recensés dans n'importe quelle surface dispo
        allParticipants = {};
        for s = 1:numel(surfaceList)
            surf = surfaceList{s};
            if isfield(S.(groupName), surf)
                Tsurf = S.(groupName).(surf);
                if ~isempty(Tsurf) && any(strcmp(Tsurf.Properties.VariableNames, 'Participant'))
                    allParticipants = [allParticipants; unique(Tsurf.Participant)]; %#ok<AGROW>
                end
            end
        end
        allParticipants = unique(allParticipants);

        for p = 1:numel(allParticipants)
            pid = allParticipants{p};

            for s = 1:numel(surfaceList)
                surf = surfaceList{s};
                if ~isfield(S.(groupName), surf); continue; end

                Tsurf = S.(groupName).(surf);
                if isempty(Tsurf); continue; end
                if ~any(strcmp(Tsurf.Properties.VariableNames, varName)); continue; end

                idx = strcmp(Tsurf.Participant, pid);
                if ~any(idx); continue; end

                y = Tsurf.(varName)(idx);
                y = y(1);

                if isempty(y) || isnan(y); continue; end

                r = r + 1;
                rows{r,1} = pid;           %#ok<AGROW>
                rows{r,2} = groupName;     %#ok<AGROW>
                rows{r,3} = surf;          %#ok<AGROW>
                rows{r,4} = y;             %#ok<AGROW>
            end
        end
    end

    if r == 0
        T_long = table();
        return;
    end

    T_long = cell2table(rows, 'VariableNames', {'Participant','Groupe','Surface','Y'});

    % Catégories ordonnées (utile pour interprétation et cohérence)
    T_long.Participant = categorical(T_long.Participant);
    T_long.Groupe  = categorical(T_long.Groupe,  groupList,  'Ordinal', true);
    T_long.Surface = categorical(T_long.Surface, surfaceList,'Ordinal', true);
end

function [p, F] = getTermStats_STD(aov, termName)
    p = NaN; F = NaN;
    if any(strcmp(aov.Term, termName))
        idx = strcmp(aov.Term, termName);
        p = aov.pValue(idx);
        F = aov.FStat(idx);
    end
end

function grid = allcombGrid(groupList, surfaceList)
% Grille complète Groupe x Surface pour EMMs
    [G,S] = ndgrid(1:numel(groupList), 1:numel(surfaceList));
    grid = table();
    grid.Groupe  = categorical(groupList(G(:))', groupList, 'Ordinal', true);
    grid.Surface = categorical(surfaceList(S(:))', surfaceList, 'Ordinal', true);

    % Participant requis dans predict() ? Non, si Conditional=false, mais
    % certains objets l'aiment présent -> on met un dummy.
    grid.Participant = categorical(repmat("DUMMY", height(grid), 1));
end

function T = pairwiseMainEffect(lme, groupList, surfaceList, effectName, alpha)
% Pairwise sur un effet principal via EMMs (population-level)
    grid = allcombGrid(groupList, surfaceList);
    mu = predict(lme, grid, 'Conditional', false);

    T = table();
    if effectName == "Groupe"
        % moyenne marginale par Groupe (moyenne sur Surface)
        groups = categorical(groupList, groupList, 'Ordinal', true);
        m = zeros(numel(groups),1);
        for i = 1:numel(groups)
            m(i) = mean(mu(grid.Groupe == groups(i)));
        end
        [T, ~] = pairwiseFromMeans(m, string(categories(groups)), "Groupe", alpha);
    elseif effectName == "Surface"
        surfs = categorical(surfaceList, surfaceList, 'Ordinal', true);
        m = zeros(numel(surfs),1);
        for i = 1:numel(surfs)
            m(i) = mean(mu(grid.Surface == surfs(i)));
        end
        [T, ~] = pairwiseFromMeans(m, string(categories(surfs)), "Surface", alpha);
    end
end

function T = pairwiseSimpleEffects(lme, groupList, surfaceList, mode, alpha)
% Simple effects via EMMs (population-level) :
% - "SurfaceWithinGroup" : compare surfaces à l'intérieur de chaque groupe
% - "GroupWithinSurface" : compare groupes à l'intérieur de chaque surface
    grid = allcombGrid(groupList, surfaceList);
    mu = predict(lme, grid, 'Conditional', false);

    T = table();

    groups = categorical(groupList, groupList, 'Ordinal', true);
    surfs  = categorical(surfaceList, surfaceList, 'Ordinal', true);

    if mode == "SurfaceWithinGroup"
        for g = 1:numel(groups)
            m = zeros(numel(surfs),1);
            for s = 1:numel(surfs)
                m(s) = mu(grid.Groupe==groups(g) & grid.Surface==surfs(s));
            end
            [Tout, ~] = pairwiseFromMeans(m, string(categories(surfs)), ...
                "Surface|Groupe=" + string(groups(g)), alpha);
            T = [T; Tout]; %#ok<AGROW>
        end

    elseif mode == "GroupWithinSurface"
        for s = 1:numel(surfs)
            m = zeros(numel(groups),1);
            for g = 1:numel(groups)
                m(g) = mu(grid.Groupe==groups(g) & grid.Surface==surfs(s));
            end
            [Tout, ~] = pairwiseFromMeans(m, string(categories(groups)), ...
                "Groupe|Surface=" + string(surfs(s)), alpha);
            T = [T; Tout]; %#ok<AGROW>
        end
    end
end

function [T, padj] = pairwiseFromMeans(m, labels, effectLabel, alpha)
% Pairwise contrasts sur des "means" m (déjà marginales / EMMs)
% Ici, on fait des différences simples + p-values via t approx :
% IMPORTANT: pour une approche strictement inférentielle, il faudrait propager
% la matrice de covariance des EMMs. MATLAB ne donne pas un "emmeans" natif,
% donc on fournit une solution "standard pratique" souvent acceptée:
% on reporte Estimate + p ajustées Holm, et on recommande de compléter par IC
% des EMMs déjà fournis (CI_Low/High) et par des plots.
%
% Si tu veux l'inférence exacte EMM (SE/DF), je peux te donner la version
% complète avec construction de la matrice de contraste à partir de lme.DesignMatrix.

    n = numel(m);
    pairs = n*(n-1)/2;

    comp1 = strings(pairs,1);
    comp2 = strings(pairs,1);
    est   = zeros(pairs,1);
    pval  = nan(pairs,1);

    k = 0;
    for i = 1:n-1
        for j = i+1:n
            k = k + 1;
            comp1(k) = labels(i);
            comp2(k) = labels(j);
            est(k)   = m(i) - m(j);

            % p-value "placeholder" = NaN (à éviter si tu veux une inférence stricte)
            % -> On laisse NaN par transparence plutôt que d'inventer un SE.
            pval(k) = NaN;
        end
    end

    % On n'invente pas des p-values sans SE/DF.
    % Mais tu veux des post-hocs "clean": la bonne solution est la version contrastes
    % sur le modèle (coefTest) via une matrice de contraste construite sur la grille EMM.
    %
    % Pour rester 100% correct, on renvoie ici uniquement les estimates.
    padj = pval;

    T = table();
    T.Effect = repmat(string(effectLabel), pairs, 1);
    T.Level1 = comp1;
    T.Level2 = comp2;
    T.Estimate = est;
    T.pValue_raw = pval;
    T.pValue_adj = padj;
    T.Signif = repmat("", pairs, 1);

    % Marqueur si p-values disponibles (sinon vide)
    if all(~isnan(padj))
        T.Signif(padj < alpha) = "*";
    end
end

function s = starsCol(p)
    s = strings(numel(p),1);

    for i = 1:numel(p)
        if isnan(p(i))
            s(i) = "";
            continue;
        end

        nStar = sum(p(i) < [0.05 0.01 0.001]);

        if nStar > 0
            s(i) = string(repmat('*', 1, nStar));  % <-- FIX ICI
        else
            s(i) = "";
        end
    end
end

function p_adj = holmAdjustVector(p)
    p = p(:);
    [ps, idx] = sort(p);
    m = numel(p);

    adj_sorted = min(ps .* (m - (1:m)' + 1), 1);

    % monotonicity
    for i = m-1:-1:1
        adj_sorted(i) = min(adj_sorted(i), adj_sorted(i+1));
    end

    p_adj = NaN(m,1);
    p_adj(idx) = adj_sorted;
end