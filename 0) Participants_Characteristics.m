%% Create_or_Update_ParticipantsMetadata.m
% Source : participants_metadata.mat (table)
% Export : participants_metadata.csv

clc; clear; close all;

% === CHEMIN DE SAUVEGARDE ===
save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\Statistical_Analysis_LMM\Stats_Descriptives';

mat_file = fullfile(save_path, 'participants_metadata.mat');
csv_file = fullfile(save_path, 'participants_metadata.csv');

% === CHARGER OU CREER ===
if isfile(mat_file)
    load(mat_file, 'T');
    fprintf('Fichier participants_metadata.mat chargé.\n');
else
    T = table( ...
    strings(0,1), categorical(strings(0,1)), nan(0,1), categorical(strings(0,1)), nan(0,1), nan(0,1), nan(0,1), nan(0,1), ...
    'VariableNames', {'Participant','AgeGroup','AgeMonths','Sex','Height_cm','Weight_kg','L0_m','IMC'} );
    fprintf('Nouvelle table participants_metadata créée.\n');
end

% === AJOUT / MISE A JOUR : PLUSIEURS PARTICIPANTS ===
newRows = table( ...
    ["CTL_01"; "CTL_02"; "CTL_03"; "CTL_04"; "CTL_05"; "CTL_06"; "CTL_07"; "CTL_08"; "CTL_09"; "CTL_10";
    "CTL_11"; "CTL_12"; "CTL_13"; "CTL_14"; "CTL_15"; "CTL_16"; "CTL_17"; "CTL_18"; "CTL_19"; "CTL_20";
    "CTL_22"; "CTL_23"; "CTL_24"; "CTL_25"; "CTL_26"; "CTL_27"; "CTL_28"; "CTL_29"; "CTL_30";
    "CTL_31"; "CTL_32"; "CTL_33"; "CTL_34"; "CTL_35"; "CTL_36"; "CTL_37"; "CTL_38"; "CTL_39"; "CTL_40";
    "CTL_41"; "CTL_42"; "CTL_44"; "CTL_46"; "CTL_47"; "CTL_48"; "CTL_49"; "CTL_50";
    "CTL_51"; "CTL_53"; "CTL_56"; "CTL_57"; "CTL_58"; "CTL_59"; "CTL_60";
    "CTL_61"; "CTL_62"; "CTL_63"; "CTL_65"; "CTL_66"; "CTL_67"; "CTL_68"; "CTL_69"; "CTL_70";
    "CTL_71"; "CTL_73"; "CTL_75"; "CTL_76"; "CTL_77";
], ...
    categorical(["Enfants"; "Enfants"; "Adultes"; "Adultes"; "Adultes"; "Enfants"; "Adultes"; "Adultes"; "Enfants"; "Enfants"; ...
    "Enfants"; "Adultes"; "Adultes"; "Jeunes Enfants"; "Jeunes Enfants"; "Enfants"; "Adolescents"; "Adultes"; "Enfants"; "Enfants"; ...
    "Adultes"; "Jeunes Enfants"; "Adultes"; "Adultes"; "Adolescents"; "Adultes"; "Adultes"; "Adultes"; "Adultes"; ...
    "Enfants"; "Adultes"; "Adultes"; "Enfants"; "Adolescents"; "Adultes"; "Jeunes Enfants"; "Enfants"; "Enfants"; "Jeunes Enfants"; ...
    "Adultes"; "Adultes"; "Enfants"; "Enfants"; "Adolescents"; "Adolescents"; "Enfants"; "Adolescents"; ...
    "Enfants"; "Jeunes Enfants"; "Adolescents"; "Enfants"; "Enfants"; "Enfants"; "Adolescents"; ...
    "Adolescents"; "Adolescents"; "Jeunes Enfants"; "Jeunes Enfants"; "Adolescents"; "Jeunes Enfants"; "Adolescents"; "Adolescents"; "Jeunes Enfants"; ...
    "Adolescents"; "Adolescents"; "Jeunes Enfants"; "Adolescents"; "Jeunes Enfants"]), ...
    ...
    [117; 87; 295; 341; 283; 96; 268; 329; 73; 134; ...
    132; 293; 264; 55; 44; 77; 169; 267; 79; 110; ...
    310; 71; 253; 350; 189; 269; 359; 226; 303; ...
    73; 295; 300; 140; 161; 338; 64; 100; 87; 58; ...
    281; 279; 98; 92; 189; 178; 127; 171; ...
    135; 62; 150; 83; 103; 103; 189; ...
    145; 161; 49; 44; 158; 53; 162; 144; 70; ...
    182; 153; 65; 181; 31;
    ], ...
    ...
    categorical(["M"; "M"; "F"; "M"; "F"; "F"; "M"; "M"; "M"; "F"; 
        "F"; "F"; "F"; "F"; "F"; "M"; "M"; "M"; "M"; "F";
        "F"; "F"; "F"; "F"; "F"; "F"; "F"; "M"; "M";
        "F"; "F"; "M"; "F"; "F"; "M"; "F"; "F"; "F"; "M";
        "F"; "M"; "F"; "F"; "F"; "F"; "M"; "F";
        "F"; "F"; "F"; "F"; "F"; "F"; "F";
        "F"; "M"; "M"; "F"; "M"; "F"; "M"; "F"; "M";
        "F"; "F"; "F"; "M"; "M";
        ]), ...
        ...
    [138; 121.5; 160; 174; 164; 135; 184; 174; 114; 150;
        150; 176; 166; 101; 85; 109; 180; 186; 113; 137;
        182; 110.5; 167; 165; 164; 167; 168; 184; 182;
        126; 176; 176; 149; 160; 181; 114; 134; 125; 109; 
        165; 188.5; 131; 130; 172; 166; 166; 178;
        154; 110; 156; 119; 123; 127; 163;
        152; 168.5; 103; 101; 164; 106; 155; 156.5; 125;
        169; 164; 113; 188; 92;
        ], ...
        ...
    [29.5; 21; 67; 78.5; 54.5; 28.7; 80; 69.1; 20.4; 43; 
        38; 72.5; 83.6; 17; 14.2; 18.5; 60.4; 77.5; 18; 30.7;
        70; 18.3; 72.1; 66.8; 60.4; 48; 54; 73; 75; 
        32; 62.9; 75.5; 44.5; 47.8; 101; 19.2; 25.7; 23.5; 19;
        56; 70.5; 30; 28; 54.6; 82; 61.3; 71; 
        44; 19.8; 47; 22; 23.2; 26.1; 50;
        35; 68; 16.8; 16; 43; 18.2; 40; 44.2; 22.5;
        53.8; 49; 20; 59.4; 15;
        ], ...
        ...
    [0.7325; 0.6675; 0.815; 0.929; 0.8825; 0.711; 0.945; 0.9025; 0.54; 0.8;
0.8125; 0.9505; 0.9; 0.475; 0.455; 0.555; 0.9725; 0.985; 0.5625; 0.7375;
0.9925; 0.55; 0.8825; 0.8735; 0.87; 0.8665; 0.9125; 0.945; 0.975;
0.6175; 0.931; 0.9; 0.793; 0.84; 0.952; 0.5645; 0.6565; 0.623; 0.54;
0.851; 1; 0.6825; 0.646; 0.8925; 0.924; 0.905; 0.9485;
0.84; 0.5385; 0.803; 0.5755; 0.608; 0.62; 0.875;
0.8045; 0.9275; 0.47; 0.472; 0.89; 0.5; 0.83; 0.8145; 0.6275;
0.8825; 0.8775; 0.568; 1.01; 0.45;
], ...
...
    'VariableNames', {'Participant','AgeGroup','AgeMonths','Sex','Height_cm','Weight_kg','L0_m'} );

% === Harmoniser les colonnes newRows
if ~ismember("IMC", string(newRows.Properties.VariableNames))
    newRows.IMC = nan(height(newRows),1);
end

% Option robuste : forcer le même ordre de colonnes que T
newRows = newRows(:, T.Properties.VariableNames);

for i = 1:height(newRows)
    pid = newRows.Participant(i);
    idx = find(T.Participant == pid, 1);

    if isempty(idx)
        T = [T; newRows(i,:)];
        fprintf('Ajout participant %s.\n', pid);
    else
        T(idx,:) = newRows(i,:);
        fprintf('Mise à jour participant %s.\n', pid);
    end
end

% === NETTOYAGE / TRI ===
T = sortrows(T, "Participant");

T.IMC = T.Weight_kg ./ ((T.Height_cm/100).^2);

% === SAUVEGARDES ===
save(mat_file, 'T');
writetable(T, csv_file);

fprintf('\nSauvegarde OK.\n');
fprintf('MAT: %s\n', mat_file);
fprintf('CSV: %s\n', csv_file);
fprintf('N participants: %d\n', height(T));