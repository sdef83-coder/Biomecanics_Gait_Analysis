%% EXTRACTION VARIABLES SPATIO-TEMPORELLES PAR GROUPE D'ÂGE ET CONDITION (format.mat et .csv) + Visualisation
% extraction radar plot plus bas : 2ème partie du script
% extraction graphique valeurs en fonction de l'âge : 3ème partie du script

clc;
clear;
close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL')
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'))
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script'))

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
            'NormWalkRatio','vitCadencePasParMinute','NormStepLength','NormCadence', ...
            'MoS_AP_Mean','MoS_ML_Mean','MoS_HS_AP','MoS_HS_ML'};

newNames = {'Single support time (%)', 'Double support time (%)','Stride width (cm)', ...
            'Gait speed (m.s^{-1})','Stride length (m)', 'Stride time (s)', ...
            'Norm WR (ua)','Cadence (step.min^{-1})', ...
            'Norm Step length (ua)', 'Norm Cadence (ua)', ...
            'MoS AP Stance (mm)', 'MoS ML Stance (mm)', 'MoS AP HS (mm)', 'MoS ML HS (mm)'};

renameMap = containers.Map(oldNames, newNames);
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
                numericCols = varfun(@isnumeric, individualData, 'output', 'uniform');
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

% Groupes d'âge (certains peuvent être absents)
groupList = {'JeunesEnfants', 'Enfants', 'Adolescents', 'Adultes'};

originalNames = oldNames;

renameMap = containers.Map(originalNames, newNames);

% Préfixes à traiter : Moyenne + CV
prefixes = {'Mean_', 'CV_'};

% CORRECTION : Structure d'export par participant avec 1 ligne par participant
for p = 1:length(prefixes)
    prefix = prefixes{p};
    
    for iVar = 1:length(originalNames)
        varTech = originalNames{iVar};
        varNameReadable = renameMap(varTech);
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

%% RADAR PLOTS 5 DOMAINES - INTER-GROUPES
clc; clear; close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));
load('SpatioTemporalDATA.mat');

save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Fig\SpatioTempo-DATA';
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

conditions = {'Plat', 'Medium', 'High'};
groups = {'JeunesEnfants', 'Enfants', 'Adolescents', 'Adultes'};

% === RADAR INTER-GROUPES (compare les groupes d'âge) ===
fprintf('\n=== Génération des radar plots INTER-GROUPES ===\n');
for i = 1:length(conditions)
    condition = conditions{i};
    fig = radarGaitPlot_5Domains_Inter(SpatioTemporalDATA, condition, groups);
    
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1400, 1000]);
    filename = fullfile(save_path, sprintf('RadarPlot_5Domains_Inter_%s.png', condition));
    print(fig, filename, '-dpng', '-r300');
    fprintf('  ✅ %s\n', condition);
    close(fig);
end

% RADAR PLOTS 5 DOMAINES - INTRA-GROUPES
fprintf('\n=== Génération des radar plots INTRA-GROUPES ===\n');

condColors = {[0.2 0.4 1], [0 0.6 0], [1 0 0]}; % Plat, Medium, High

for g = 1:length(groups)
    groupName = groups{g};
    fig = radarGaitPlot_5Domains_Intra(SpatioTemporalDATA, groupName, conditions, condColors);
    
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1400, 1000]);
    filename = fullfile(save_path, sprintf('RadarPlot_5Domains_Intra_%s.png', groupName));
    print(fig, filename, '-dpng', '-r300');
    fprintf('  ✅ %s\n', groupName);
    close(fig);
end

fprintf('\n✅ TOTAL: 7 radar plots générés (3 inter + 4 intra)\n');
fprintf('📂 Sauvegardés dans: %s\n', save_path);

%% NUAGE DE POINTS DE L'EVOLUTION DES PARAMETRES SPATIO-TEMPORELLES EN FONCTION DU TEMPS
clc, clear, close all;

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result')
load('SpatioTemporalDATA.mat')

% Créer le dossier de sortie si nécessaire
output_folder = fullfile('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Fig', 'SpatioTempo-DATA');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% === DÉFINITION DES VARIABLES À TRACER (organisées par domaine) ===
variables_to_plot = {
    % PACE
    'Mean_Gait speed (m.s^{-1})', 'Mean_Stride length (m)', 'Mean_Stride time (s)', ...
    'CV_Gait speed (m.s^{-1})', 'CV_Stride length (m)', 'CV_Stride time (s)', ...
    
    % RHYTHM
    'Mean_Cadence (step.min^{-1})', 'Mean_Normalized Cadence (ua)', 'Mean_Normalized Walk ratio (ua)', ...
    
    % STABILITY
    'Mean_Stride width (cm)', 'Mean_MoS HS ML (mm)', 'Mean_Double support time (%)', ...
    'Mean_Single support time (%)', ...
    'Mean_MoS AP Mean (mm)', 'Mean_MoS ML Mean (mm)', 'Mean_MoS HS AP (mm)', ...
    
    % ASYMMETRY
    'SI_Stride time (s)', 'SI_Stride length (m)', 'SI_Stride width (cm)', ...
    
    % VARIABILITY
    'CV_Stride time (s)', 'CV_Stride length (m)', 'CV_Stride width (cm)', ...
    'CV_Double support time (%)', 'CV_Single support time (%)'
};

% Couleurs pour les 3 surfaces
color_map = containers.Map(...
    {'Plat', 'Medium', 'High'}, ...
    {[0 0.447 0.741], [0 0.6 0], [0.85 0.1 0.1]}); % Bleu, Vert, Rouge

% Fusionner toutes les conditions dans une seule table
DATA_all = [SpatioTemporalDATA.ALL.Plat; ...
            SpatioTemporalDATA.ALL.Medium; ...
            SpatioTemporalDATA.ALL.High];

% Définir les tranches d'âge (en mois)
tranches = [36, 72; 72, 144; 144, 216; 216, 432];
nTranches = size(tranches, 1);
tranche_names = {'Young Children', 'Children', 'Adolescents', 'Adults'};

% Boucle sur chaque variable
for i = 1:length(variables_to_plot)
    varname = variables_to_plot{i};
    
    % Vérifier si la variable existe dans les données
    if ~any(strcmp(DATA_all.Properties.VariableNames, varname))
        fprintf('⚠️  Variable non trouvée : %s\n', varname);
        continue;
    end

    fig = figure('Position', [100, 100, 1200, 700]);
    hold on;

    % Lignes verticales pour les tranches d'âge
    for t = 1:nTranches
        xline(tranches(t,1), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'HandleVisibility','off');
    end
    xline(tranches(end,2), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'HandleVisibility','off');

    % Structure pour stocker les points d'inflexion
    inflection_info = struct();
    
    % Boucle sur les surfaces
    conditions = {'Plat', 'Medium', 'High'};
    for c = 1:length(conditions)
        cond_name = conditions{c};
        color = color_map(cond_name);
        data_cond = DATA_all(strcmp(DATA_all.Condition, cond_name), :);

        % --- Nuage de points individuels
        scatter(data_cond.AgeMonths, data_cond.(varname), 20, ...
            'filled', ...
            'MarkerFaceColor', color, ...
            'MarkerEdgeColor', color, ...
            'MarkerFaceAlpha', 0.25, ...
            'MarkerEdgeAlpha', 0.25, ...
            'DisplayName', cond_name);

        % --- Moyennes et SD par tranche
        moyennes = nan(nTranches,1);
        SD = nan(nTranches,1);
        x_center = nan(nTranches,1);

        for iT = 1:nTranches
            borne_inf = tranches(iT,1);
            borne_sup = tranches(iT,2);
            idx = data_cond.AgeMonths >= borne_inf & data_cond.AgeMonths < borne_sup;
            x_vals = data_cond.AgeMonths(idx);
            y_vals = data_cond.(varname)(idx);

            if ~isempty(y_vals)
                moyennes(iT) = mean(y_vals, 'omitnan');
                SD(iT) = std(y_vals, 'omitnan');
                x_center(iT) = mean(x_vals, 'omitnan');
            end
        end

        % --- Tracer moyennes avec barres d'erreur
        errorbar(x_center, moyennes, SD, '-', ...
            'Color', color, 'LineWidth', 2, ...
            'Marker', 'o', 'MarkerFaceColor', color, ...
            'CapSize', 6, 'DisplayName', ['Mean ' cond_name]);

        % --- SPLINE LISSEE SUR LES POINTS INDIVIDUELS (plus robuste)
        x = data_cond.AgeMonths;
        y = data_cond.(varname);
        valid = ~isnan(x) & ~isnan(y);
        x = x(valid); y = y(valid);
        [x, ord] = sort(x); y = y(ord);

        try
            if numel(x) >= 10
                p_smooth = 0.90;                 % 0.85–0.95 selon le bruit
                sp = csaps(x, y, p_smooth);      % spline de lissage

                % Courbe lissée pour visuel
                xx = linspace(min(x), max(x), 300);
                yy = fnval(sp, xx);
                plot(xx, yy, '-', 'Color', min(color*0.8 + 0.2, 1), 'LineWidth', 2.5, ...
                    'DisplayName', [cond_name ' (spline)']);

                % === DÉTECTION DU POINT D'INFLEXION (2e dérivée = 0)
                sp2 = fnder(sp, 2);
                z = fnzeros(sp2, [min(x) max(x)]);

                % Normaliser la sortie de fnzeros en liste de candidats
                if isempty(z)
                    x0 = [];
                elseif isvector(z)
                    x0 = z(:).';
                else
                    x0 = mean(z,1); % centres des intervalles de zéro
                end

                % Filtre: éviter les bords (artefacts fréquents)
                if ~isempty(x0)
                    margin = 0.05*(max(x)-min(x));
                    x0 = x0(x0 > (min(x)+margin) & x0 < (max(x)-margin));
                end

                if ~isempty(x0)
                    % Intensité d'inflexion via 3e dérivée
                    sp3 = fnder(sp, 3);
                    [~, ixbest] = max(abs(fnval(sp3, x0)));

                    x_inf = x0(ixbest);
                    y_inf = fnval(sp, x_inf);

                    % Marquer le point d'inflexion
                    plot(x_inf, y_inf, 'p', ...
                        'MarkerSize', 12, 'MarkerFaceColor', color, ...
                        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
                        'HandleVisibility', 'off');

                    % Rattacher à la tranche correspondante pour ton annotation
                    t_idx = find(tranches(:,1) <= x_inf & x_inf < tranches(:,2), 1, 'first');
                    if isempty(t_idx)
                        % si x_inf est à la toute fin, le rattacher à la dernière tranche
                        t_idx = nTranches;
                    end

                    inflection_info.(cond_name).method       = 'spline-2nd-deriv';
                    inflection_info.(cond_name).x_pos        = x_inf;
                    inflection_info.(cond_name).value        = y_inf;
                    inflection_info.(cond_name).tranche      = t_idx;
                    inflection_info.(cond_name).age_range    = tranches(t_idx, :);
                    inflection_info.(cond_name).tranche_name = tranche_names{t_idx};
                end
            end
        catch ME
            % Fallback propre si Spline Toolbox absent ou autre souci
            warning('Spline/derivative step failed for %s (%s): %s', cond_name, varname, ME.message);
        end

    end % fin boucle conditions

    % === ANNOTATION DES POINTS D'INFLEXION ===
    if ~isempty(fieldnames(inflection_info))
        annotation_text = {'\bfInflection Points:'};
        y_offset = 0.88; % Position verticale de départ
        
        for c = 1:length(conditions)
            cond_name = conditions{c};
            if isfield(inflection_info, cond_name)
                info = inflection_info.(cond_name);
                color = color_map(cond_name);

                text_line = sprintf('\\color[rgb]{%.2f %.2f %.2f}%s: %s (%.0f-%.0f mo) — x\\_inf=%.0f mo', ...
                    color(1), color(2), color(3), ...
                    cond_name, info.tranche_name, ...
                    info.age_range(1), info.age_range(2), ...
                    info.x_pos);

                annotation_text{end+1} = text_line; %#ok<AGROW>
            end
        end
        
        annotation('textbox', [0.62 y_offset 0.35 0.12], ...
            'String', annotation_text, ...
            'FitBoxToText', 'on', ...
            'BackgroundColor', 'white', ...
            'EdgeColor', [0.3 0.3 0.3], ...
            'LineWidth', 1, ...
            'FontSize', 9, ...
            'Interpreter', 'tex');
    end

    % === MISE EN FORME ===
    xlabel('Age (months)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(strrep(varname, '_', ' '), 'Interpreter', 'tex', 'FontSize', 12, 'FontWeight', 'bold');
    title(strrep(varname, '_', ' '), 'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    legend('Location', 'eastoutside', 'Box', 'off', 'FontSize', 10);
    grid on; box on;
    ax = gca; ax.GridAlpha = 0.3; ax.GridLineStyle = ':';

    % Sauvegarde
    safe_varname = regexprep(varname, '[^\w]', '_');
    saveas(gcf, fullfile(output_folder, [safe_varname '_vs_Age_Inflection.png']));
    close(gcf);
    
    fprintf('✅ Figure générée : %s\n', varname);
end

disp("✅ Toutes les figures générées avec points d'inflexion détectés.");