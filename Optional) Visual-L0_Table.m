clc, clear, close all;

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL')

load("l0_participants.mat", 'l0_map')

% Récupère les clés et les valeurs
keys_list = keys(l0_map);
values_list = values(l0_map);

% Crée une table avec les noms et les longueurs
T = table(keys_list', cell2mat(values_list)', ...
    'VariableNames', {'Participant', 'L0_m'});

% Affiche le tableau
disp(T);

