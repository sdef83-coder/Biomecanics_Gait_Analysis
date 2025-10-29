%% EXTRACTION VARIABLES SPATIO-TEMPORELLES PAR GROUPE D'ÂGE ET CONDITION (format.mat et .csv) + Visualisation
% extraction radar plot plus bas : 2ème partie du script
% extraction graphique valeurs en fonction de l'âge : 3ème partie du script

clc;
clear;
close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL')
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'))
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script'))

% Dossier où sont sauvegardés les résultats MoS par participant
mos_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\MoS';

% Chemin de sauvegarde
save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result';
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% Attribution des participants dans des groupes
ParticipantGroup;
Condition = {'Plat', 'Medium', 'High'};
% === Dictionnaire de renommage des variables ===
oldNames = {'pctSimpleAppuie', 'DoubleSupport','LargeurPas','vitFoulee','distFoulee', 'tempsFoulee', ...
            'NormWalkRatio','vitCadencePasParMinute','NormStepLength','NormCadence'};

newNames = {'Single support time (%)', 'Double support time (%)','Stride width (cm)', ...
            'Gait speed (m.s^{-1})','Stride length (m)', 'Stride time (s)', ...
            'Norm WR (ua)','Cadence (step.min^{-1})', ...
            'Norm Step length (ua)', 'Norm Cadence (ua)'};

renameMap = containers.Map(oldNames, newNames);

% Ajoute après la définition de renameMap
renameMap('MoS_AP_HS_mm')   = 'MoS AP HS (mm)';
renameMap('MoS_AP_mean_mm') = 'MoS AP Stance (mm)';
renameMap('MoS_ML_HS_mm')   = 'MoS ML HS (mm)';
renameMap('MoS_ML_mean_mm') = 'MoS ML Stance (mm)';

% si tu exportes aussi %L0 :
renameMap('MoS_AP_HS_pL0')   = 'MoS AP HS (%L0)';
renameMap('MoS_AP_mean_pL0') = 'MoS AP Stance (%L0)';
renameMap('MoS_ML_HS_pL0')   = 'MoS ML HS (%L0)';
renameMap('MoS_ML_mean_pL0') = 'MoS ML Stance (%L0)';

prefixes = {'Mean_', 'CV_'};

% Initialisation
SpatioTemporalDATA = struct();
groupNames = fieldnames(Group);

% Première phase : collecte des données individuelles (incluant médianes)
for g = 1:length(groupNames)
    groupName = groupNames{g};
    participants = Group.(groupName);
    fprintf('Traitement du groupe : %s\n', groupName);
    
    for iC = 1:length(Condition)
        cond = Condition{iC};
        fprintf('  Condition : %s\n', cond);
        
        % Préallocation d’un tableau de structs
        % recapRows = repmat(struct(), 1, length(participants));
        % rowIdx = 1;
        rowIdx = 1;
        recapRows = [];

        for iP = 1:length(participants)
            participant = participants{iP};
            fprintf('    Participant traité : %s\n', participant);
            file = [participant '_' cond '.mat'];

            if exist(file, 'file')
                data = load(file);
                stats = Spatiotempocalc(data.c.resultsAll.kin.Left, data.c.resultsAll.kin.Right);
                DATA.(participant).(cond).stats = stats;

                row = struct();
                row.Participant = string(participant);
                row.Condition = string(cond);
                row.AgeMonths = Association_Age(participant);

                % ==== Récup MoS pour ce participant et cette condition ====
try
    mosAgg = get_mos_aggregates(participant, cond, mos_dir);
    
    % Bruts (mm)
    row.Mean_MoS_AP_HS_mm   = mosAgg.MoS_AP_HS_mm;
    row.Mean_MoS_AP_mean_mm = mosAgg.MoS_AP_mean_mm;
    row.Mean_MoS_ML_HS_mm   = mosAgg.MoS_ML_HS_mm;
    row.Mean_MoS_ML_mean_mm = mosAgg.MoS_ML_mean_mm;

    % %L0 si dispo dans les fichiers
    if isfield(mosAgg,'MoS_AP_HS_pL0')
        row.Mean_MoS_AP_HS_pL0   = mosAgg.MoS_AP_HS_pL0;
        row.Mean_MoS_AP_mean_pL0 = mosAgg.MoS_AP_mean_pL0;
        row.Mean_MoS_ML_HS_pL0   = mosAgg.MoS_ML_HS_pL0;
        row.Mean_MoS_ML_mean_pL0 = mosAgg.MoS_ML_mean_pL0;
    end

    % Garde aussi dans DATA (optionnel)
    DATA.(participant).(cond).MoS = mosAgg;

catch ME
    warning('MoS manquant pour %s - %s : %s', participant, cond, ME.message);
end
% ==== fin Récup MoS ===

                % Nombre de cycles
                nCyclesLeft = size(data.c.resultsAll.kin.Left, 2);
                nCyclesRight = size(data.c.resultsAll.kin.Right, 2);
                row.NCycles_Left = nCyclesLeft;
                row.NCycles_Right = nCyclesRight;

                DATA.(participant).(cond).nCycles.Left = nCyclesLeft;
                DATA.(participant).(cond).nCycles.Right = nCyclesRight;

                statsFields = fieldnames(stats);
                for f = 1:length(statsFields)
                    fname = statsFields{f};

                    if endsWith(fname, '_Mean_Mean')
                        shortName = extractBefore(fname, '_Mean_Mean');
                        row.(['Mean_' shortName]) = stats.(fname);
                    elseif endsWith(fname, '_CV_Mean')
                        shortName = extractBefore(fname, '_CV_Mean');
                        row.(['CV_' shortName]) = stats.(fname);
                    elseif endsWith(fname, '_SI')
                        shortName = extractBefore(fname, '_SI');
                        row.(['SI_' shortName]) = stats.(fname);
                    end
                end

                   % Initialisation correcte du tableau struct si premier passage
        if rowIdx == 1
            recapRows = repmat(row, 1, length(participants));
        end

        recapRows(rowIdx) = row;
        rowIdx = rowIdx + 1;

if ~isfield(SpatioTemporalDATA, 'ALL') || ~isfield(SpatioTemporalDATA.ALL, cond)
    SpatioTemporalDATA.ALL.(cond) = row; 
else
    SpatioTemporalDATA.ALL.(cond)(end+1) = row;
end

            else
                warning('Fichier manquant : %s', file);
            end
        end

        % Convertir en table une fois à la fin de la condition
        recapData = struct2table(recapRows(1:rowIdx-1));
        SpatioTemporalDATA.(groupName).(cond) = recapData;
        fprintf('  -> Données du groupe %s pour la condition %s traitées.\n', groupName, cond);
    end
end


% Deuxième phase : Créer les matrices de comparaison par condition (incluant médianes)
disp('--- Début de la génération des matrices de comparaison par condition ---');
for iC = 1:length(Condition)
    cond = Condition{iC};
    
    % Collecter toutes les variables numériques communes
    allVariables = {};
    
    % Identifier les variables à partir du premier groupe qui a des données
    for g = 1:length(groupNames)
        groupName = groupNames{g};
        if isfield(SpatioTemporalDATA, groupName) && isfield(SpatioTemporalDATA.(groupName), cond)
            if ~isempty(SpatioTemporalDATA.(groupName).(cond))
                individualData = SpatioTemporalDATA.(groupName).(cond);
                numericCols = varfun(@isnumeric, individualData, 'OutputFormat', 'uniform');
                allVariables = individualData.Properties.VariableNames(numericCols);
                break;
            end
        end
    end
    
    if ~isempty(allVariables)
        % Créer les noms de colonnes : Mean_Var1, Median_Var1, Std_Var1, Mean_Var2, Median_Var2, Std_Var2...
        columnNames = {};
        for v = 1:length(allVariables)
            varName = allVariables{v};
            columnNames{end+1} = ['Mean_' varName];
            columnNames{end+1} = ['Std_' varName];
            columnNames{end+1} = ['Median_' varName];
        end
        
        % Initialiser la matrice de comparaison
        comparisonMatrix = array2table(NaN(length(groupNames), length(columnNames)), ...
            'VariableNames', columnNames, ...
            'RowNames', groupNames);
        
        % Remplir la matrice pour chaque groupe
        for g = 1:length(groupNames)
            groupName = groupNames{g};
            
            if isfield(SpatioTemporalDATA, groupName) && isfield(SpatioTemporalDATA.(groupName), cond)
                if ~isempty(SpatioTemporalDATA.(groupName).(cond))
                    individualData = SpatioTemporalDATA.(groupName).(cond);
                    
                    % Calculer statistiques pour chaque variable
                    for v = 1:length(allVariables)
                        varName = allVariables{v};
                        
                        if ismember(varName, individualData.Properties.VariableNames)
                            values = individualData.(varName);
                            cleanValues = values(~isnan(values));
                            
                            if ~isempty(cleanValues)
                                comparisonMatrix{groupName, ['Mean_' varName]} = mean(cleanValues);
                                comparisonMatrix{groupName, ['Std_' varName]} = std(cleanValues);
                                comparisonMatrix{groupName, ['Median_' varName]} = median(cleanValues);
                            end
                        end
                    end
                end
            end
        end
        
        % Sauvegarder la matrice de comparaison au niveau supérieur
        SpatioTemporalDATA.(cond) = comparisonMatrix;
        fprintf('  → Matrice de comparaison générée pour la condition : %s\n', cond);
    end
end

for iC = 1:length(Condition)
    cond = Condition{iC};
    if isfield(SpatioTemporalDATA.ALL, cond)
        SpatioTemporalDATA.ALL.(cond) = struct2table(SpatioTemporalDATA.ALL.(cond));
    end
end

% Export des matrices de comparaison (incluant médianes)
comparisonFile = fullfile(save_path, 'Comparaison_Groupes_SpatioTemporel.xlsx');
for iC = 1:length(Condition)
    cond = Condition{iC};
    if isfield(SpatioTemporalDATA, cond)
        comparisonMatrix = SpatioTemporalDATA.(cond);
        writetable(comparisonMatrix, comparisonFile, 'Sheet', cond, 'WriteRowNames', true);
    end
end

% === RENOMMAGE DES VARIABLES POUR LES TABLEAUX INDIVIDUELS (ALL) ===
for iC = 1:length(Condition)
    cond = Condition{iC};
    if isfield(SpatioTemporalDATA.ALL, cond)
        tableData = SpatioTemporalDATA.ALL.(cond);
        vars = tableData.Properties.VariableNames;

        for iV = 1:length(vars)
            varName = vars{iV};

            prefixes = {'Mean_', 'CV_', 'SI_'};
            matchedPrefix = '';
            baseName = varName;

            for p = 1:length(prefixes)
                if startsWith(varName, prefixes{p})
                    matchedPrefix = prefixes{p};
                    baseName = extractAfter(varName, matchedPrefix);
                    break;
                end
            end

            if isKey(renameMap, baseName)
                newBase = renameMap(baseName);
                tableData.Properties.VariableNames{iV} = [matchedPrefix newBase];
            end
        end

        SpatioTemporalDATA.ALL.(cond) = tableData;
    end
end

% === RENOMMAGE DES VARIABLES POUR CHAQUE GROUPE ===
groupList = {'JeunesEnfants', 'Enfants', 'Adolescents', 'Adultes'}; % adapte selon tes groupes présents

for iC = 1:length(Condition)
    cond = Condition{iC};
    
    for g = 1:length(groupList)
        groupName = groupList{g};
        
        if isfield(SpatioTemporalDATA, groupName) && isfield(SpatioTemporalDATA.(groupName), cond)
            tableData = SpatioTemporalDATA.(groupName).(cond);
            vars = tableData.Properties.VariableNames;

            for iV = 1:length(vars)
                varName = vars{iV};

                prefixes = {'Mean_', 'CV_', 'SI_'};
                matchedPrefix = '';
                baseName = varName;

                for p = 1:length(prefixes)
                    if startsWith(varName, prefixes{p})
                        matchedPrefix = prefixes{p};
                        baseName = extractAfter(varName, matchedPrefix);
                        break;
                    end
                end

                if isKey(renameMap, baseName)
                    newBase = renameMap(baseName);
                    tableData.Properties.VariableNames{iV} = [matchedPrefix newBase];
                end
            end

            % Remettre dans la structure renommée
            SpatioTemporalDATA.(groupName).(cond) = tableData;
        end
    end
end

% === EXPORT INDIVIDUEL DES VARIABLES PAR CSV ===
% Dossier de sortie
csv_export_path = fullfile(save_path, 'Matrice', 'CSV_Variables');
if ~exist(csv_export_path, 'dir')
    mkdir(csv_export_path);
end

% Groupes d'âge 
groupList = {'JeunesEnfants', 'Enfants', 'Adolescents', 'Adultes'};

% Base spatio-temporelle
originalNames = oldNames;
newNamesBase  = newNames;  % même ordre que oldNames

% Ajoute les variables MoS à exporter (si tu veux aussi %L0, garde ces lignes)
mosTech = {'MoS_AP_HS_mm','MoS_AP_mean_mm','MoS_ML_HS_mm','MoS_ML_mean_mm', ...
           'MoS_AP_HS_pL0','MoS_AP_mean_pL0','MoS_ML_HS_pL0','MoS_ML_mean_pL0'};
mosReadable = {'MoS AP HS (mm)','MoS AP Stance (mm)','MoS ML HS (mm)','MoS ML Stance (mm)', ...
               'MoS AP HS (%L0)','MoS AP Stance (%L0)','MoS ML HS (%L0)','MoS ML Stance (%L0)'};

% Fusionne proprement
originalNames = [originalNames, mosTech];
newNamesAll   = [newNamesBase, mosReadable];

% Map OK (même longueur des deux côtés)
renameMapExport = containers.Map(originalNames, newNamesAll);

% Préfixes à traiter : Moyenne + CV
prefixes = {'Mean_', 'CV_'};

% CORRECTION : Structure d'export par participant avec 1 ligne par participant
for p = 1:length(prefixes)
    prefix = prefixes{p};
    
    for iVar = 1:length(originalNames)
        varTech = originalNames{iVar};
        varNameReadable = renameMapExport(varTech);
        varFull = [prefix varNameReadable];  % ex: Mean_Gait speed (m.s^{-1})

        % Initialisation - CORRECTION ICI
        exportData = table();
        exportRowIndex = 1;
        
        % Obtenir tous les participants uniques de tous les groupes
        allParticipants = {};
        participantGroups = {};
        
        % Première passe : identifier tous les participants
        for g = 1:length(groupList)
            gName = groupList{g};
            
            % Vérifier si le groupe existe et a des données
            if isfield(SpatioTemporalDATA, gName)
                % Prendre n'importe quelle condition pour obtenir la liste des participants
                for iC = 1:length(Condition)
                    cond = Condition{iC};
                    if isfield(SpatioTemporalDATA.(gName), cond)
                        T = SpatioTemporalDATA.(gName).(cond);
                        if ~isempty(T) && any(strcmp(T.Properties.VariableNames, varFull))
                            groupParticipants = unique(T.Participant);
                            for i = 1:length(groupParticipants)
                                participantID = groupParticipants{i};
                                if ~ismember(participantID, allParticipants)
                                    allParticipants{end+1} = participantID;
                                    participantGroups{end+1} = gName;
                                end
                            end
                        end
                        break; % On n'a besoin que d'une condition pour avoir la liste
                    end
                end
            end
        end
        
        % Deuxième passe : collecter les données pour chaque participant
        for i = 1:length(allParticipants)
            participantID = allParticipants{i};
            gName = participantGroups{i};
            
            % Initialiser la ligne pour ce participant
            newRow = table();
            newRow.GroupeAge = {gName};
            newRow.Participant = {participantID};
            newRow.Plat = NaN;
            newRow.Medium = NaN;
            newRow.High = NaN;
            
            % Collecter les valeurs pour chaque condition
            for iC = 1:length(Condition)
                cond = Condition{iC};
                
                if isfield(SpatioTemporalDATA, gName) && isfield(SpatioTemporalDATA.(gName), cond)
                    T = SpatioTemporalDATA.(gName).(cond);
                    
                    if ~isempty(T) && any(strcmp(T.Properties.VariableNames, varFull))
                        % Trouver la ligne correspondant à ce participant
                        participantRows = strcmp(T.Participant, participantID);
                        
                        if any(participantRows)
                            value = T.(varFull)(participantRows);
                            if ~isempty(value) && ~isnan(value(1))
                                newRow.(cond) = value(1);
                            end
                        end
                    end
                end
            end
            
            % Ajouter la ligne au tableau d'export
            if exportRowIndex == 1
                exportData = newRow;
            else
                exportData = [exportData; newRow];
            end
            exportRowIndex = exportRowIndex + 1;
        end
        
        % Sauvegarder le fichier CSV
        if ~isempty(exportData)
            % Nettoyer le nom du fichier
            safeFileName = matlab.lang.makeValidName([prefix '_' varNameReadable]);
            csvPath = fullfile(csv_export_path, [safeFileName '.csv']);
            writetable(exportData, csvPath);
            fprintf('✅ Export format large : %s (%d participants)\n', csvPath, height(exportData));
        else
            fprintf('⚠️  Aucune donnée trouvée pour : %s\n', varFull);
        end
    end
end

% OPTIONNEL : Afficher un résumé des exports
fprintf('\n=== RÉSUMÉ DES EXPORTS ===\n');
fprintf('Format: 1 ligne par participant avec colonnes Plat, Medium, High\n');
for g = 1:length(groupList)
    gName = groupList{g};
    if isfield(SpatioTemporalDATA, gName)
        % Compter les participants de ce groupe
        participantCount = 0;
        for iC = 1:length(Condition)
            cond = Condition{iC};
            if isfield(SpatioTemporalDATA.(gName), cond)
                T = SpatioTemporalDATA.(gName).(cond);
                if ~isempty(T)
                    participantCount = length(unique(T.Participant));
                    break;
                end
            end
        end
        fprintf('%s: %d participants\n', gName, participantCount);
    end
end

% Sauvegarder la structure finale
save(fullfile(save_path, 'SpatioTemporalDATA.mat'), 'SpatioTemporalDATA');

% Export des tableaux individuels (ALL) en format CSV pour chaque surface
for iC = 1:length(Condition)
    cond = Condition{iC};
    if isfield(SpatioTemporalDATA.ALL, cond)
        tableToExport = SpatioTemporalDATA.ALL.(cond);
        % Nom du fichier CSV
        csvFileName = fullfile(save_path, ['SpatioTemporal_ALL_' cond '.csv']);
        % Export en CSV
        writetable(tableToExport, csvFileName);
        fprintf('→ Export CSV réalisé pour : %s\n', cond);
    end
end

% Affichage des résultats
disp('=== RÉSULTATS ===');
disp('Structure SpatioTemporalDATA sauvegardée avec succès.');

%% Fonction utilitaire
function mosAgg = get_mos_aggregates(participant, cond, mos_dir)
% Retourne des moyennes par condition pour les indicateurs MoS clés
% - suppose un fichier: mos_dir/MoS_results_<participant>.mat
% - lit la table MoS_data.results et filtre Surface==cond

    f = fullfile(mos_dir, sprintf('MoS_results_%s.mat', participant));
    if ~exist(f,'file')
        error('Fichier MoS introuvable: %s', f);
    end

    S = load(f, 'MoS_data');
    if ~isfield(S,'MoS_data') || ~istable(S.MoS_data.results)
        error('Structure MoS_data.results invalide pour %s', participant);
    end
    T = S.MoS_data.results;

    % Filtre par condition
    if ~ismember('Surface', T.Properties.VariableNames)
        error('La table MoS ne contient pas la colonne Surface.');
    end
    idx = strcmp(T.Surface, cond);
    T = T(idx, :);
    if isempty(T)
        error('Aucun cycle pour la condition %s.', cond);
    end

    % Champs requis (brut mm)
    req = {'MoS_Heel_Strike_AP','MoS_AP_Mean','MoS_Heel_Strike_ML','MoS_ML_Mean'};
    for k = 1:numel(req)
        if ~ismember(req{k}, T.Properties.VariableNames)
            error('Champ manquant dans MoS: %s', req{k});
        end
    end

    % Moyennes brutes (mm)
    mosAgg.MoS_AP_HS_mm   = mean(T.MoS_Heel_Strike_AP, 'omitnan');
    mosAgg.MoS_AP_mean_mm = mean(T.MoS_AP_Mean,        'omitnan');
    mosAgg.MoS_ML_HS_mm   = mean(T.MoS_Heel_Strike_ML, 'omitnan');
    mosAgg.MoS_ML_mean_mm = mean(T.MoS_ML_Mean,        'omitnan');

    % Versions %L0 si elles existent
    optP = {'MoS_Heel_Strike_AP_P','MoS_AP_Mean_P','MoS_Heel_Strike_ML_P','MoS_ML_Mean_P'};
    haveP = all(ismember(optP, T.Properties.VariableNames));
    if haveP
        mosAgg.MoS_AP_HS_pL0   = mean(T.MoS_Heel_Strike_AP_P, 'omitnan');
        mosAgg.MoS_AP_mean_pL0 = mean(T.MoS_AP_Mean_P,        'omitnan');
        mosAgg.MoS_ML_HS_pL0   = mean(T.MoS_Heel_Strike_ML_P, 'omitnan');
        mosAgg.MoS_ML_mean_pL0 = mean(T.MoS_ML_Mean_P,        'omitnan');
    end
end