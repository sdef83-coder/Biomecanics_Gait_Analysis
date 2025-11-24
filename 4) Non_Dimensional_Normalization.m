%% AJOUT DE PARAMETRE DANS LE .MAT POUR SPATIO-TEMPORAL ANALYSIS
% Pas obligé de runner tout les participants

clc;
clear;
close all;

% Configuration générale
conditions = {'Plat', 'Medium', 'High'};
folder_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL';
cd(folder_path);

% Chargement des valeurs L0 pour tous les participants
load('l0_participants.mat', 'l0_map');
g = 9.81;

% Liste des participants à traiter (peut être automatiquement générée)
%participants = keys(l0_map); % Récupère automatiquement tous les participants du fichier L0
participants = {'CTL_65', 'CTL_31', 'CTL_38', 'CTL_39', 'CTL_44', 'CTL_46', 'CTL_03', 'CTL_13', 'CTL_22', 'CTL_29'}; % Si besoin de changer qu'un seul participant, décommentez et commentez la ligne au dessus

% Boucle pour chaque participant
for p = 1:length(participants)
    participant = participants{p};
    
    % Récupération de la valeur L0 spécifique au participant
    if isKey(l0_map, participant)
        l0 = l0_map(participant);
        fprintf('Traitement du participant %s avec L0 = %.4f\n', participant, l0);
    else
        fprintf('Participant %s non trouvé dans les données L0, passage au suivant\n', participant);
        continue;
    end
    
    % Boucle pour chaque condition
    for i = 1:length(conditions)
        condition = conditions{i};
        filename = [participant '_' condition '.mat'];
        
        if isfile(filename)
            load(filename); % charge la structure'c'
            
            % === Jambe droite : ajout NormStepLength, NormCadence, NormWalkRatio, NormStepWidthHeel, et NormWalkSpeed ===
            N_right = length(c.resultsAll.kin.Right);
            for n = 1:N_right
                distPasR = c.resultsAll.kin.Right(n).distPas;
                cadenceR = c.resultsAll.kin.Right(n).vitCadencePasParMinute * 2; % x2 car ici cycle par minute et pas ppm
                LargPasR = c.resultsAll.kin.Right(n).stepWidthHeel * 10; % x10 car ici la variable est en mm, on l'a met en cm
                WalkSpeedR = c.resultsAll.kin.Right(n).vitFoulee;
                % (Hof et al., 1996) pour les méthodes de normalisation des variables
                c.resultsAll.kin.Right(n).NormStepLength = distPasR / l0;
                c.resultsAll.kin.Right(n).NormCadence = (cadenceR/60) * sqrt(l0 / g);
                c.resultsAll.kin.Right(n).NormWalkRatio = (distPasR / l0) / ((cadenceR/60) * sqrt(l0 / g));
                c.resultsAll.kin.Right(n).NormStepWidthHeel = LargPasR / l0;
                c.resultsAll.kin.Right(n).NormWalkSpeed = WalkSpeedR / sqrt(l0*g);
            end
            
            % === Jambe gauche : ajout NormStepLength, NormCadence, NormWalkRatio, NormStepWidthHeel, et NormWalkSpeed ===
            N_left = length(c.resultsAll.kin.Left);
            for n = 1:N_left
                distPasL = c.resultsAll.kin.Left(n).distPas;
                cadenceL = c.resultsAll.kin.Left(n).vitCadencePasParMinute * 2;
                LargPasL = c.resultsAll.kin.Left(n).stepWidthHeel * 10;
                WalkSpeedL = c.resultsAll.kin.Left(n).vitFoulee;
                % (Hof et al., 1996) pour les méthodes de normalisation des variables
                c.resultsAll.kin.Left(n).NormStepLength = distPasL / l0;
                c.resultsAll.kin.Left(n).NormCadence = (cadenceL/60) * sqrt(l0 / g);
                c.resultsAll.kin.Left(n).NormWalkRatio = (distPasL / l0) / ((cadenceL/60) * sqrt(l0 / g));
                c.resultsAll.kin.Left(n).NormStepWidthHeel = LargPasL / l0;
                c.resultsAll.kin.Left(n).NormWalkSpeed = WalkSpeedL / sqrt(l0*g);
            end
            
            % === Sauvegarde du fichier avec les données ajoutées ===
            save(filename, 'c', '-v7.3');
            fprintf('  → Données ajoutées à : %s\n', filename);
        else
            fprintf('  → Fichier non trouvé : %s\n', filename);
        end
    end
    fprintf('Participant %s terminé\n\n', participant);
end

fprintf('Traitement terminé pour tous les participants!\n');