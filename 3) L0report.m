%% Renseigner le L0 de l'ensemble des participants pour la suite du Script
% L0 = Longueur de jambe moyenne entre gauche et droite (en mètre)
clc, clear, close all;

% Chemin de sauvegarde
save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL';

% Charger la Map existante ou créer une nouvelle
full_filename = fullfile(save_path, 'l0_participants.mat');
if isfile(full_filename)
    load(full_filename, 'l0_map');
    fprintf('Fichier L0 existant chargé\n');
else
    l0_map = containers.Map();
    fprintf('Nouvelle Map L0 créée\n');
end

% === PARTICIPANTS EXISTANTS (à commenter après la première exécution) ===
% Décommentez SEULEMENT si je dois recréer le fichier complet ou modifier des valeurs existantes
% l0_map('CTL_01') = 0.7325;
% l0_map('CTL_02') = 0.6675;
% l0_map('CTL_03') = 0.815 ;
% l0_map('CTL_04') = 0.929;
% l0_map('CTL_05') = 0.8825;
% l0_map('CTL_06') = 0.711;
% l0_map('CTL_07') = 0.945;
% l0_map('CTL_08') = 0.9025;
% l0_map('CTL_09') = 0.54;
% l0_map('CTL_10') = 0.80;
% l0_map('CTL_11') = 0.8125;
% l0_map('CTL_12') = 0.9505;
% l0_map('CTL_13') = 0.90 ;
% l0_map('CTL_14') = 0.475;
% l0_map('CTL_15') = 0.455;
% l0_map('CTL_16') = 0.555;
% l0_map('CTL_17') = 0.9725;
% l0_map('CTL_18') = 0.985;
% l0_map('CTL_19') = 0.5625;
% l0_map('CTL_20') = 0.7375;
% l0_map('CTL_22') = 0.9925;
% l0_map('CTL_23') = 0.55;
% l0_map('CTL_24') = 0.8825;
% l0_map('CTL_25') = 0.8735;
% l0_map('CTL_26') = 0.87;
% l0_map('CTL_27') = 0.8665;
% l0_map('CTL_28') = 0.9125;
% l0_map('CTL_29') = 0.945;
% l0_map('CTL_30') = 0.975;
% l0_map('CTL_31') = 0.6175;
% l0_map('CTL_32') = 0.931;
% l0_map('CTL_33') = 0.900;
% l0_map('CTL_34') = 0.793;
% l0_map('CTL_35') = 0.84;
% l0_map('CTL_36') = 0.952;
% l0_map('CTL_37') = 0.5645;
% l0_map('CTL_38') = 0.6565;
% l0_map('CTL_39') = 0.623;
% l0_map('CTL_40') = 0.54;
% l0_map('CTL_41') = 0.851;
% l0_map('CTL_42') = 1.00;
% l0_map('CTL_44') = 0.6825;
% l0_map('CTL_46') = 0.646;
% l0_map('CTL_47') = 0.8925;
% l0_map('CTL_48') = 0.924;
% l0_map('CTL_49') = 0.905;
% l0_map('CTL_50') = 0.9485 ;
% l0_map('CTL_51') = 0.840 ;
% l0_map('CTL_53') = 0.5385 ;
% l0_map('CTL_56') = 0.803 ;
% l0_map('CTL_57') = 0.5755 ;
% l0_map('CTL_58') = 0.608 ;
% l0_map('CTL_59') = 0.620 ;
% l0_map('CTL_60') = 0.875 ;
% l0_map('CTL_61') = 0.8045 ;
% l0_map('CTL_62') = 0.9275 ;
% l0_map('CTL_63') = 0.47;
% l0_map('CTL_65') = 0.472;
% l0_map('CTL_66') = 0.89 ;
% l0_map('CTL_67') = 0.50;
% l0_map('CTL_68') = 0.83 ;
% l0_map('CTL_69') = 0.8145;
% l0_map('CTL_70') = 0.6275;
% l0_map('CTL_71') = 0.8825;
% l0_map('CTL_73') = 0.8775;
% l0_map('CTL_75') = 0.568 ;
% l0_map('CTL_76') = 1.01;
% l0_map('CTL_77') = 0.45;

% === NOUVEAUX PARTICIPANTS (à ajouter) ===

l0_map('CTL_77') = 0.45; % Nouveau participant - décommentez si besoin de rajouter un participant


% l0_map('CTL_21') = 0.90 ; Si besoin
% l0_map('CTL_43') = 0.92 ;
% l0_map('CTL_45') = 0.9625 ;

% Sauvegarder dans le répertoire spécifié
save(full_filename, 'l0_map');

% Afficher un résumé
participants_list = keys(l0_map);
fprintf('Fichier l0_participants.mat mis à jour avec succès!\n');
fprintf('Chemin: %s\n', save_path);
fprintf('Nombre total de participants: %d\n', length(participants_list));
fprintf('Participants inclus: %s\n', strjoin(sort(participants_list), ', '));
fprintf('\n💡 PROCHAINE ÉTAPE : Lancer Non_Dimensional_Normalization.m après avoir vérifié si participant(s) dans ParticipantGroup.m\n');