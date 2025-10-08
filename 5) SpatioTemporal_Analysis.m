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
oldNames = {'pctSimpleAppuie','DoubleSupport','LargeurPas','vitFoulee','distFoulee',...
            'NormWalkRatio','vitCadencePasParMinute','NormStepLength','NormCadence'};

newNames = {'Single support time (%)','Double support time (%)','Stride width (mm)', ...
            'Gait speed (m.s^{-1})','Stride length (m)', ...
            'Normalized Walk ratio (ua)','Cadence (step.min^{-1})', ...
            'Normalized Step length (ua)', 'Normalized Cadence (ua)'};

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

% Noms techniques d'origine
originalNames = {'pctSimpleAppuie','DoubleSupport','LargeurPas','vitFoulee','distFoulee',...
                 'NormWalkRatio','vitCadencePasParMinute','NormStepLength','NormCadence'};

% Noms lisibles pour l'export
newNames = {'Single support time (%)','Double support time (%)','Stride width (mm)', ...
            'Gait speed (m.s^{-1})','Stride length (m)', ...
            'Normalized Walk ratio (ua)','Cadence (step.min^{-1})', ...
            'Normalized Step length (ua)', 'Normalized Cadence (ua)'};

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

%% RADAR PLOT DES PARAMETRES SPATIOTEMPORELS A LA MARCHE
clc; clear; close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));
load('SpatioTemporalDATA.mat');

% Dossier de sauvegarde
save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Fig\SpatioTempo-DATA';

% Créer le dossier s'il n'existe pas
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% Paramètres étudiés
param_names = {'Gait speed (m.s^{-1})','Stride length (m)', 'Normalized Step length (ua)', ...
    'Normalized Walk ratio (ua)','Cadence (step.min^{-1})', 'Normalized Cadence (ua)', ...
            'Single support time (%)','Double support time (%)', 'Stride width (mm)'};
conditions = {'Plat', 'Medium', 'High'};

for i = 1:length(conditions)
    condition = conditions{i};
    fig = radarGaitPlot(SpatioTemporalDATA, condition, param_names);
    
    % Définir la taille de la figure avant sauvegarde
    set(fig, 'Units', 'pixels', 'Position', [100, 100, 1200, 800]);
    
    % Sauvegarder la figure avec haute résolution
    filename = fullfile(save_path, sprintf('RadarPlot_%s.png', condition));
    print(fig, filename, '-dpng', '-r300');
    fprintf('Figure sauvegardée : %s\n', filename);
end

%%% === INTRAGROUPE JEUNES ENFANTS === %%%
group = 'JeunesEnfants';
colors = {[0.2 0.4 1], [0 0.6 0], [1 0 0]}; % bleu, vert, rouge

fig_jeunes_mean = makeIntragroupRadar(SpatioTemporalDATA, group, conditions, param_names, colors);
set(fig_jeunes_mean, 'Units', 'pixels', 'Position', [100, 100, 1200, 800]);
filename = fullfile(save_path, 'RadarPlot_Intragroupe_JeunesEnfants.png');
print(fig_jeunes_mean, filename, '-dpng', '-r300');
fprintf('Figure sauvegardée : %s\n', filename);

%%% === INTRAGROUPE ENFANTS === %%%
group = 'Enfants';

fig_enfants_mean = makeIntragroupRadar(SpatioTemporalDATA, group, conditions, param_names, colors);
set(fig_enfants_mean, 'Units', 'pixels', 'Position', [100, 100, 1200, 800]);
filename = fullfile(save_path, 'RadarPlot_Intragroupe_Enfants.png');
print(fig_enfants_mean, filename, '-dpng', '-r300');
fprintf('Figure sauvegardée : %s\n', filename);

%%% === INTRAGROUPE ADOLESCENTS === %%%
group = 'Adolescents';

fig_adolescents_mean = makeIntragroupRadar(SpatioTemporalDATA, group, conditions, param_names, colors);
set(fig_adolescents_mean, 'Units', 'pixels', 'Position', [100, 100, 1200, 800]);
filename = fullfile(save_path, 'RadarPlot_Intragroupe_Adolescents.png');
print(fig_adolescents_mean, filename, '-dpng', '-r300');
fprintf('Figure sauvegardée : %s\n', filename);

%%% === INTRAGROUPE ADULTES === %%%
group = 'Adultes';

fig_adultes_mean = makeIntragroupRadar(SpatioTemporalDATA, group, conditions, param_names, colors);
set(fig_adultes_mean, 'Units', 'pixels', 'Position', [100, 100, 1200, 800]);
filename = fullfile(save_path, 'RadarPlot_Intragroupe_Adultes.png');
print(fig_adultes_mean, filename, '-dpng', '-r300');
fprintf('Figure sauvegardée : %s\n', filename);

fprintf('\nToutes les figures ont été sauvegardées dans :\n%s\n', save_path);

%% NUAGE DE POINTS DE L'EVOLUTION DES PARAMETRES SPATIO-TEMPORELLES EN FONCTION DU TEMPS
clc, clear, close all;

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result')
load('SpatioTemporalDATA.mat')

% Créer le dossier de sortie si nécessaire
output_folder = fullfile('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Fig', 'SpatioTempo-DATA');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Variables à tracer
variables_to_plot = {
    'Mean_Single support time (%)','Mean_Double support time (%)','Mean_Stride width (mm)', ...
    'Mean_Gait speed (m.s^{-1})','Mean_Stride length (m)', ...
    'Mean_Normalized Walk ratio (ua)','Mean_Cadence (step.min^{-1})', ...
    'Mean_Normalized Step length (ua)', 'Mean_Normalized Cadence (ua)', ...
    'CV_Single support time (%)','CV_Double support time (%)','CV_Stride width (mm)', ...
    'CV_Gait speed (m.s^{-1})','CV_Stride length (m)', ...
    'CV_Normalized Walk ratio (ua)','CV_Cadence (step.min^{-1})', ...
    'CV_Normalized Step length (ua)', 'CV_Normalized Cadence (ua)'};

% Couleurs pour les 3 surfaces
color_map = containers.Map(...
    {'Plat', 'Medium', 'High'}, ...
    {[0 0.447 0.741], [0 0.6 0], [0.85 0.1 0.1]}); % Bleu, Vert, Rouge

% Fusionner toutes les conditions dans une seule table
DATA_all = [SpatioTemporalDATA.ALL.Plat; ...
            SpatioTemporalDATA.ALL.Medium; ...
            SpatioTemporalDATA.ALL.High];

% Définir les tranches d’âge (en mois)
tranches = [36, 72; 72, 144; 144, 216; 216, 432];
nTranches = size(tranches, 1);

% Boucle sur chaque variable
for i = 1:length(variables_to_plot)
    varname = variables_to_plot{i};

    figure;
    hold on;

    % Lignes verticales pour les tranches d’âge
    for t = 1:nTranches
        xline(tranches(t,1), '--k', 'LineWidth', 1.2, 'HandleVisibility','off');
    end
    xline(tranches(end,2), '--k', 'LineWidth', 1.2, 'HandleVisibility','off');

    % Boucle sur les surfaces
    for cond = {'Plat', 'Medium', 'High'}
        cond_name = cond{1};
        color = color_map(cond_name);
        data_cond = DATA_all(strcmp(DATA_all.Condition, cond_name), :);

        % Tracer les points individuels
        scatter(data_cond.AgeMonths, data_cond.(varname), 20, ...
    'filled', ...
    'MarkerFaceColor', color, ...
    'MarkerEdgeColor', color, ...
    'MarkerFaceAlpha', 0.25, ...
    'MarkerEdgeAlpha', 0.25, ...
    'DisplayName', cond_name);

        % Moyennes et SD par tranche
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

        % Tracer moyennes avec barres d’erreur
        errorbar(x_center, moyennes, SD, '-', ...
            'Color', color, 'LineWidth', 2, ...
            'Marker', 'o', 'MarkerFaceColor', color, ...
            'CapSize', 6, 'DisplayName', ['Moyenne ' cond_name]);
    end

    % Mise en forme
    xlabel('Âge (mois)', 'FontSize', 12);
    ylabel(varname, 'Interpreter','none', 'FontSize', 12);
   title({[strrep(varname, '_', ' ')]}, 'FontSize', 13);
    legend('Location','eastoutside', 'Box', 'off');
    grid on;
    box on;

    % Sauvegarde
    saveas(gcf, fullfile(output_folder, [regexprep(varname, '[^\w]', '_') '_vs_Age_Moyennes.png']));
    close;
end

disp("✅ Figures générées avec moyennes ± SD par tranche et par surface.");