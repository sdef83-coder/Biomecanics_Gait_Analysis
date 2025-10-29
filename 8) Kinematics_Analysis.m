%% EXTRACTION DES MATRICES SUR LA CINEMATIQUE DES MEMBRES INF. SELON LES 3 PLANS DE L'ESPACE

clc;
clear;
close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL')
addpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI');
% Définir le chemin pour l'enregistrement des résultats
save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Fig';
% Vérifier si le dossier existe, sinon le créer
if ~exist(save_path, 'dir')
    mkdir(save_path);
    disp(['Création du dossier de sauvegarde: ' save_path]);
end
% Définir le chemin d'enregistrement des CSV
csv_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Matrice\CSV';
% Charger les groupes de participants
groupe_a_etudier = 'Adolescents'; % 'Enfants', 'Adolescents', 'Adultes',
ParticipantGroup;
Participant = Group.(groupe_a_etudier);
Condition = {'Plat' 'Medium' 'High'};
Plan = {'Sagittal', 'Frontal', 'Transverse'};

for iP = 1:length(Participant)
    for iC = 1:length(Condition)
        file = [Participant{iP} '_' Condition{iC} '.mat'];
        data = load(file);
        
        % Extraction des données de chaque articulation dans les trois plans
        % Plan sagittal (colonne 1), frontal (colonne 2) et transverse (colonne 3)
        for iPlane = 1:length(Plan)
            plane = Plan{iPlane};
            
            % Extraction en utilisant le format correct: data.c.results.MeanLeg.angAtFullCycle.Ankle{1, 1}(:,2)
            % Vérifier que les colonnes existent avant de les extraire
            try
                % Extraire les données pour le plan spécifique (numéro de colonne correspondant)
                DATA.(Participant{iP}).(Condition{iC}).(plane).Ankle = data.c.results.MeanLeg.angAtFullCycle.Ankle{1, 1}(:, iPlane);
                DATA.(Participant{iP}).(Condition{iC}).(plane).Knee = data.c.results.MeanLeg.angAtFullCycle.Knee{1, 1}(:, iPlane);
                DATA.(Participant{iP}).(Condition{iC}).(plane).Hip = data.c.results.MeanLeg.angAtFullCycle.Hip{1, 1}(:, iPlane);
                
                % Afficher des informations diagnostiques (seulement pour le premier fichier)
                if iP == 1 && iC == 1 && iPlane == 1
                    disp(['Chargement réussi des données pour le fichier: ' file]);
                    disp(['Taille des données Ankle pour le plan ' plane ': ' num2str(size(DATA.(Participant{iP}).(Condition{iC}).(plane).Ankle))]);
                end
            catch
                % En cas d'erreur, utiliser des NaN et afficher un message d'avertissement
                % Déterminer la taille des données en utilisant le premier plan disponible si possible
                if iPlane == 1
                    warning(['Impossible d''accéder aux données pour ' Participant{iP} ', ' Condition{iC} ', ' plane]);
                    continue;
                else
                    try
                        data_size = size(data.c.results.MeanLeg.angAtFullCycle.Ankle{1, 1}(:, 1));
                        DATA.(Participant{iP}).(Condition{iC}).(plane).Ankle = nan(data_size);
                        DATA.(Participant{iP}).(Condition{iC}).(plane).Knee = nan(data_size);
                        DATA.(Participant{iP}).(Condition{iC}).(plane).Hip = nan(data_size);
                    catch
                        warning(['Structure de données inattendue pour ' Participant{iP} ', ' Condition{iC}]);
                        continue;
                    end
                end
            end
        end
    end
end

% Définition des articulations à traiter
Articulations = {'Ankle', 'Knee', 'Hip'};

% Définition des couleurs pour chaque condition
couleurs = {'b', 'g', 'r'}; % Bleu pour Plat, Vert pour Medium, Rouge pour High

% Création d'un vecteur pour l'axe X (0-100% du cycle de marche)
x = linspace(0, 100, size(DATA.(Participant{1}).(Condition{1}).Sagittal.Ankle, 1));

% Création des matrices pour chaque articulation et chaque condition dans chaque plan
% où chaque ligne représente un sujet et chaque colonne un point du cycle
disp('Création des matrices par articulation, condition et plan...');

% Déterminer la taille des données pour initialiser les matrices
nb_points = size(DATA.(Participant{1}).(Condition{1}).Sagittal.Ankle, 1);

% Initialiser les matrices pour toutes les combinaisons articulation/condition/plan
for iPlane = 1:length(Plan)
    plane = Plan{iPlane};
    
    for iA = 1:length(Articulations)
        articulation = Articulations{iA};
        
        for iC = 1:length(Condition)
            condition = Condition{iC};
            
            % Créer une matrice où chaque ligne est un sujet et chaque colonne un point du cycle
            % Initialiser avec NaN pour marquer les données manquantes
            MATRICES.(plane).(articulation).(condition) = nan(length(Participant), nb_points);
            
            % Remplir la matrice avec les données de chaque participant
            for iP = 1:length(Participant)
                participant = Participant{iP};
                
                % MODIFICATION: Exclure la cheville de CTL_19 des calculs
                if strcmp(participant, 'CTL_19') && strcmp(articulation, 'Ankle')
                    % Ne pas inclure les données de la cheville pour CTL_19
                    % Laisser la ligne avec NaN pour exclure ce participant de la moyenne
                    disp(['Exclusion de la cheville pour ' participant ' dans le plan ' plane ' et condition ' condition]);
                    continue;
                end
                
                % Vérifier si les données existent pour ce participant/condition/plan
                if isfield(DATA, participant) && isfield(DATA.(participant), condition) && ...
                   isfield(DATA.(participant).(condition), plane) && isfield(DATA.(participant).(condition).(plane), articulation)
                    
                    % Extraire les données et les transposer pour les mettre en ligne
                    participant_data = DATA.(participant).(condition).(plane).(articulation)(:, 1)';
                    MATRICES.(plane).(articulation).(condition)(iP, :) = participant_data;
                end
            end
            
            % Calculer le cycle moyen pour cette combinaison plan/articulation/condition
            % en ignorant les valeurs NaN (incluant maintenant CTL_19 pour la cheville)
            valid_data = MATRICES.(plane).(articulation).(condition);
            valid_data = valid_data(~isnan(valid_data(:,1)), :);  % Supprimer les lignes avec NaN
            
            if ~isempty(valid_data)
                MATRICES.(plane).(articulation).(['Mean_' condition]) = mean(valid_data, 1);
                MATRICES.(plane).(articulation).(['Std_' condition]) = std(valid_data, 0, 1);
                
                % Affichage du nombre de participants inclus dans le calcul
                if strcmp(articulation, 'Ankle')
                    disp(['Matrice pour ' plane ' - ' articulation ' - ' condition ': ' ...
                          num2str(size(valid_data, 1)) ' sujets inclus dans la moyenne (CTL_19 exclu) × ' ...
                          num2str(size(MATRICES.(plane).(articulation).(condition), 2)) ' points']);
                else
                    disp(['Matrice pour ' plane ' - ' articulation ' - ' condition ': ' ...
                          num2str(size(valid_data, 1)) ' sujets × ' ...
                          num2str(size(MATRICES.(plane).(articulation).(condition), 2)) ' points']);
                end
            else
                MATRICES.(plane).(articulation).(['Mean_' condition]) = nan(1, nb_points);
                MATRICES.(plane).(articulation).(['Std_' condition]) = nan(1, nb_points);
                disp(['Aucune donnée valide pour ' plane ' - ' articulation ' - ' condition]);
            end
            
            disp(['Cycle moyen pour ' plane ' - ' articulation ' - ' condition ' calculé']);
        end
    end
end

% Sauvegarder la structure MATRICES dans un fichier .mat
matrices_file = fullfile(save_path, ['MATRICES' groupe_a_etudier '.mat']);
save(matrices_file, 'MATRICES');
disp(['Structure MATRICES sauvegardée dans le fichier: ' matrices_file]);

% Définir le chemin d'enregistrement des CSV
% Définir le chemin d'enregistrement des CSV dans le sous-dossier du groupe
csv_base_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Matrice\CSV';
csv_path = fullfile(csv_base_path, groupe_a_etudier);

% Créer le dossier s'il n'existe pas
if ~exist(csv_path, 'dir')
    mkdir(csv_path);
end

% Créer le dossier s'il n'existe pas
if ~exist(csv_path, 'dir')
    mkdir(csv_path);
end

% Exporter uniquement les matrices de données individuelles avec noms de participants
for iPlane = 1:length(Plan)
    plane = Plan{iPlane};
    
    for iA = 1:length(Articulations)
        articulation = Articulations{iA};
        
        for iC = 1:length(Condition)
            condition = Condition{iC};
            
            % Extraire la matrice des données individuelles
            data_matrix = MATRICES.(plane).(articulation).(condition);
            
            % Vérifier que la taille correspond au nombre de participants
            if size(data_matrix, 1) ~= length(Participant)
                warning(['Taille incohérente pour ' articulation ' - ' plane ' - ' condition]);
                continue;
            end

            % Créer une table avec noms des participants
            T = array2table(data_matrix);
            T = addvars(T, Participant', 'Before', 1, 'NewVariableNames', 'Participant');
            
            % Créer des noms de colonnes pour les % du cycle
            nPoints = size(data_matrix, 2);
            labels = arrayfun(@(x) sprintf('P%.1f', x), linspace(0, 100, nPoints), 'UniformOutput', false);
            T.Properties.VariableNames(2:end) = labels;

            % Définir le nom de fichier CSV
            filename_csv = fullfile(csv_path, ...
                ['Individual_' articulation '_' plane '_' condition '.csv']);
            
            % Sauvegarder la table en CSV
            writetable(T, filename_csv);
        end
    end
end

disp(['✅ Export CSV dans : ' csv_path]);

% Pour chaque plan et articulation, calcul et affichage du cycle moyen et écart-type
figure_count = 0;
for iPlane = 1:length(Plan)
    plane = Plan{iPlane};
    
    for iA = 1:length(Articulations)
        articulation = Articulations{iA};
        figure_count = figure_count + 1;
        
        % Création de la figure
        fig = figure('Name', ['Cycle moyen - ' articulation ' - ' plane], 'Position', [100, 100, 800, 500]);
        hold on;
        
        % Pour chaque condition
        for iC = 1:length(Condition)
            condition = Condition{iC};
            
            % Utiliser directement les cycles moyens calculés et stockés dans la matrice
            mean_data = MATRICES.(plane).(articulation).(['Mean_' condition])';
            std_data = MATRICES.(plane).(articulation).(['Std_' condition])';

            % Tracé des courbes individuelles
individual_data = MATRICES.(plane).(articulation).(condition);
for iP = 1:size(individual_data, 1)
    % Vérifie que les données sont valides (pas NaN)
    if ~any(isnan(individual_data(iP, :)))
        plot(x, individual_data(iP, :), 'Color', couleurs{iC}, ...
            'LineWidth', 0.5, 'LineStyle', '--', 'HandleVisibility', 'off');
    end
end
            % Tracé de la moyenne
            plot(x, mean_data, 'Color', couleurs{iC}, 'LineWidth', 2, 'DisplayName', condition);
            
            % Tracé de l'écart-type (zone ombrée) sans l'ajouter à la légende
            fill([x, fliplr(x)], [mean_data + std_data; flipud(mean_data - std_data)], couleurs{iC}, ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        
        % Personnalisation du graphique
        xlabel('Cycle de marche (%)');
        ylabel('Angle (°)');
        
        % Personnalisation du titre
        title(['Cycles moyens ' groupe_a_etudier ' - ' articulation ' - Plan ' plane]);
        
        grid on;
        legend('Location', 'best');
        hold off;
        
        % Sauvegarde de la figure en format .fig et .png
        fig_filename = fullfile(save_path, ...
    [groupe_a_etudier '_' articulation '_' plane]);
        print(fig, [fig_filename '.png'], '-dpng', '-r300');
        disp(['Figure sauvegardée: ' fig_filename '.png']);
    end
end

% Affichage des informations sur le nombre de participants
disp(['Nombre de participants: ' num2str(length(Participant))]);
disp(['Toutes les figures et matrices ont été sauvegardées dans: ' save_path]);
