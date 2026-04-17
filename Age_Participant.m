%% CALCUL L'AGE DES PARTICIPANT (Années et Mois) LORS DE LA MANIP

clc, clear, close all
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Functions'));

dateNaissance = '17-02-2024';
dateEvaluation = '10-04-2026';

[annees, mois, ageEnMois] = calculAge_anneeMois(dateNaissance, dateEvaluation);

fprintf("Âge au moment de l'évaluation : %d ans et %d mois\n", annees, mois);
fprintf("Soit un total de %d mois\n", ageEnMois);


