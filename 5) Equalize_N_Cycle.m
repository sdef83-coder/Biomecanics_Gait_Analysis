%% ÉGALISATION ÉQUILIBRÉE DES CYCLES ENTRE CONDITIONS
% Pour chaque participant :
% - Même nombre TOTAL de cycles pour Plat, Medium, High
% - Répartition équilibrée Left/Right (50/50 si possible)
% - Maximise le nombre de cycles conservés
% - Sélection aléatoire (avec seed fixe pour être repro)

clc; clear; close all;

% Configuration
conditions = {'Plat', 'Medium', 'High'};
folder_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL';
cd(folder_path);

% Charger L0
load('l0_participants.mat', 'l0_map');
participants = keys(l0_map);

% === OPTION : traiter 1, plusieurs ou tous les participants ===
RUN_ONLY_SOME = true;   % false = traiter tous les participants
PARTICIPANTS_TO_RUN = {'CTL_30', 'CTL_33'};  % <--- liste de 1 ou plusieurs participant(s)

if RUN_ONLY_SOME
    % Vérifier que tous les IDs existent
    for k = 1:numel(PARTICIPANTS_TO_RUN)
        if ~ismember(PARTICIPANTS_TO_RUN{k}, participants)
            error('Le participant "%s" n''existe pas dans l0_map.', PARTICIPANTS_TO_RUN{k});
        end
    end
    % Filtrer la liste
    participants = PARTICIPANTS_TO_RUN;
    fprintf('⚙️  Mode "participants sélectionnés" activé :\n');
    disp(participants);
    fprintf('\n');
else
    fprintf('⚙️  Mode "tous les participants" activé (%d participants)\n\n', numel(participants));
end

% === PARAMÈTRES ===
RANDOM_SEED = 42;  % Seed fixe pour reproductibilité
SELECTION_METHOD = 'random'; % 'random' avec seed fixe
TARGET_LEFT_RIGHT_RATIO = 0.5; % 50% Left, 50% Right (idéal)
TOLERANCE = 0.1; % Tolérance : accepter 40%-60% pour chaque jambe

fprintf('=== ÉGALISATION ÉQUILIBRÉE DES CYCLES ===\n');
fprintf('Seed aléatoire : %d (reproductible)\n', RANDOM_SEED);
fprintf('Ratio cible Left/Right : %.0f%% / %.0f%%\n\n', ...
    TARGET_LEFT_RIGHT_RATIO*100, (1-TARGET_LEFT_RIGHT_RATIO)*100);

% Fixer le seed global
rng(RANDOM_SEED);

% Structure pour statistiques
cycle_stats = struct();

%% === TRAITEMENT PAR PARTICIPANT ===
for p = 1:length(participants)
    participant = participants{p};
    fprintf('🔄 Traitement : %s\n', participant);
    
    % === ÉTAPE 1 : Compter les cycles disponibles pour toutes les conditions ===
    cycle_counts = struct();
    files_exist = true;
    
    for iC = 1:length(conditions)
        cond = conditions{iC};
        filename = [participant '_' cond '.mat'];
        
        if ~isfile(filename)
            fprintf('  ⚠️  Fichier manquant : %s\n', filename);
            files_exist = false;
            break;
        end
        
        load(filename, 'c');
        
        cycle_counts.(cond).Left = length(c.resultsAll.kin.Left);
        cycle_counts.(cond).Right = length(c.resultsAll.kin.Right);
        cycle_counts.(cond).Total = cycle_counts.(cond).Left + cycle_counts.(cond).Right;
        
        fprintf('  📊 %s : Left=%d, Right=%d, Total=%d\n', ...
            cond, cycle_counts.(cond).Left, cycle_counts.(cond).Right, ...
            cycle_counts.(cond).Total);
    end
    
    if ~files_exist
        continue;
    end
    
    % === ÉTAPE 2 : Trouver le nombre TOTAL minimum entre les 3 conditions ===
    total_cycles = [cycle_counts.Plat.Total, cycle_counts.Medium.Total, cycle_counts.High.Total];
    min_total = min(total_cycles);
    
    fprintf('  📊 Total minimum entre conditions : %d cycles\n', min_total);
    
    % === ÉTAPE 3 : Calculer la répartition optimale Left/Right ===
    % Objectif : se rapprocher de 50/50 tout en respectant les contraintes
    
    % Calculer les disponibilités minimales par jambe
    min_left_available = min([cycle_counts.Plat.Left, cycle_counts.Medium.Left, cycle_counts.High.Left]);
    min_right_available = min([cycle_counts.Plat.Right, cycle_counts.Medium.Right, cycle_counts.High.Right]);
    
    % Stratégie : maximiser le total en équilibrant Left/Right
    [n_left_target, n_right_target] = optimize_left_right_split(...
        min_total, min_left_available, min_right_available, ...
        TARGET_LEFT_RIGHT_RATIO, TOLERANCE);
    
    fprintf('  ✂️  Cycles à conserver : Left=%d, Right=%d (Total=%d)\n', ...
        n_left_target, n_right_target, n_left_target + n_right_target);
    fprintf('  📐 Ratio obtenu : %.1f%% Left / %.1f%% Right\n', ...
        (n_left_target/(n_left_target+n_right_target))*100, ...
        (n_right_target/(n_left_target+n_right_target))*100);
    
    % Sauvegarder les stats
    cycle_stats.(participant).original = cycle_counts;
    cycle_stats.(participant).selected.Left = n_left_target;
    cycle_stats.(participant).selected.Right = n_right_target;
    cycle_stats.(participant).selected.Total = n_left_target + n_right_target;
    cycle_stats.(participant).selected.LeftRatio = n_left_target / (n_left_target + n_right_target);
    
    % === ÉTAPE 4 : Sélectionner et sauvegarder les cycles pour chaque condition ===
    for iC = 1:length(conditions)
        cond = conditions{iC};
        filename = [participant '_' cond '.mat'];
        
        load(filename, 'c');
        
        % Backup de l'original
        c.resultsAll.kin.Left_original = c.resultsAll.kin.Left;
        c.resultsAll.kin.Right_original = c.resultsAll.kin.Right;
        
        % Sélection aléatoire avec seed spécifique pour cette condition
        % (pour garantir la reproductibilité)
        participant_seed = RANDOM_SEED + p*1000 + iC*10; % Seed unique par participant/condition
        
        % Sélectionner Left
        if cycle_counts.(cond).Left > n_left_target
            selected_idx_left = select_cycles_random(...
                c.resultsAll.kin.Left, n_left_target, participant_seed);
            c.resultsAll.kin.Left = c.resultsAll.kin.Left(selected_idx_left);
            fprintf('  ✅ %s Left : %d → %d cycles (seed=%d)\n', ...
                cond, cycle_counts.(cond).Left, n_left_target, participant_seed);
        elseif cycle_counts.(cond).Left < n_left_target
            fprintf('  ⚠️  %s Left : %d cycles disponibles (< %d demandés)\n', ...
                cond, cycle_counts.(cond).Left, n_left_target);
        end
        
        % Sélectionner Right
        if cycle_counts.(cond).Right > n_right_target
            selected_idx_right = select_cycles_random(...
                c.resultsAll.kin.Right, n_right_target, participant_seed + 1);
            c.resultsAll.kin.Right = c.resultsAll.kin.Right(selected_idx_right);
            fprintf('  ✅ %s Right : %d → %d cycles (seed=%d)\n', ...
                cond, cycle_counts.(cond).Right, n_right_target, participant_seed + 1);
        elseif cycle_counts.(cond).Right < n_right_target
            fprintf('  ⚠️  %s Right : %d cycles disponibles (< %d demandés)\n', ...
                cond, cycle_counts.(cond).Right, n_right_target);
        end
        
        % Métadonnées
        c.meta.cycles_equalized = true;
        c.meta.equalization_method = 'random_balanced';
        c.meta.equalization_date = datetime('now');
        c.meta.random_seed = RANDOM_SEED;
        c.meta.participant_seed = participant_seed;
        c.meta.original_cycle_count.Left = cycle_counts.(cond).Left;
        c.meta.original_cycle_count.Right = cycle_counts.(cond).Right;
        c.meta.target_cycle_count.Left = n_left_target;
        c.meta.target_cycle_count.Right = n_right_target;
        c.meta.target_total = n_left_target + n_right_target;
        c.meta.target_left_ratio = TARGET_LEFT_RIGHT_RATIO;
        c.meta.actual_left_ratio = n_left_target / (n_left_target + n_right_target);
        
        % Sauvegarder
        save(filename, 'c', '-v7.3');
    end
    
    fprintf('  ✅ %s : Égalisation terminée\n\n', participant);
end

%% === RÉSUMÉ ET STATISTIQUES ===
fprintf('\n=== RÉSUMÉ DE L EGALISATION ===\n');

% Sauvegarder les statistiques
save('cycle_equalization_stats_balanced.mat', 'cycle_stats', 'RANDOM_SEED', ...
    'TARGET_LEFT_RIGHT_RATIO', 'TOLERANCE');

% Créer un tableau récapitulatif
summary_table = create_summary_table_balanced(cycle_stats, participants, conditions);
writetable(summary_table, 'cycle_equalization_summary_balanced.csv');

fprintf('📊 Statistiques sauvegardées dans cycle_equalization_stats_balanced.mat\n');
fprintf('📊 Tableau récapitulatif exporté dans cycle_equalization_summary_balanced.csv\n');

% Statistiques globales
fprintf('\n📈 STATISTIQUES GLOBALES :\n');
all_totals = [];
all_left_ratios = [];

for p = 1:length(participants)
    participant = participants{p};
    if isfield(cycle_stats, participant)
        all_totals = [all_totals; cycle_stats.(participant).selected.Total];
        all_left_ratios = [all_left_ratios; cycle_stats.(participant).selected.LeftRatio];
    end
end

fprintf('  Cycles totaux conservés : min=%d, max=%d, moyenne=%.1f\n', ...
    min(all_totals), max(all_totals), mean(all_totals));
fprintf('  Ratio Left moyen : %.1f%% (écart-type : %.1f%%)\n', ...
    mean(all_left_ratios)*100, std(all_left_ratios)*100);
fprintf('  Ratio cible était : %.0f%%\n', TARGET_LEFT_RIGHT_RATIO*100);

% Vérifier si les ratios respectent la tolérance
within_tolerance = abs(all_left_ratios - TARGET_LEFT_RIGHT_RATIO) <= TOLERANCE;
fprintf('  Participants dans la tolérance (±%.0f%%) : %d/%d (%.1f%%)\n', ...
    TOLERANCE*100, sum(within_tolerance), length(within_tolerance), ...
    (sum(within_tolerance)/length(within_tolerance))*100);

fprintf('\n✅ Égalisation terminée pour %d participants\n', length(participants));
fprintf('🔄 Reproductibilité garantie avec seed=%d\n', RANDOM_SEED);
fprintf('\n💡 PROCHAINE ÉTAPE : Lancer MOS.m\n');

%% ========== FONCTIONS ==========

function [n_left, n_right] = optimize_left_right_split(total_target, ...
    max_left, max_right, target_ratio, tolerance)
% Optimise la répartition Left/Right pour maximiser le nombre de cycles
% tout en respectant le ratio cible et les contraintes
%
% Entrées :
%   total_target  : nombre total de cycles à atteindre
%   max_left      : nombre maximum de cycles Left disponibles
%   max_right     : nombre maximum de cycles Right disponibles
%   target_ratio  : ratio cible pour Left (0.5 = 50%)
%   tolerance     : tolérance acceptable (0.1 = ±10%)
%
% Sorties :
%   n_left, n_right : nombre de cycles à conserver pour chaque jambe

    % Calcul idéal
    n_left_ideal = round(total_target * target_ratio);
    n_right_ideal = total_target - n_left_ideal;
    
    % Vérifier si la répartition idéale est possible
    if n_left_ideal <= max_left && n_right_ideal <= max_right
        n_left = n_left_ideal;
        n_right = n_right_ideal;
        return;
    end
    
    % Sinon, ajuster en fonction des contraintes
    if n_left_ideal > max_left
        % Pas assez de cycles Left disponibles
        n_left = max_left;
        n_right = min(total_target - n_left, max_right);
    elseif n_right_ideal > max_right
        % Pas assez de cycles Right disponibles
        n_right = max_right;
        n_left = min(total_target - n_right, max_left);
    end
    
    % Vérifier si on respecte la tolérance
    actual_ratio = n_left / (n_left + n_right);
    
    if abs(actual_ratio - target_ratio) > tolerance
        % Si hors tolérance, essayer de réajuster
        % Stratégie : réduire le total pour respecter le ratio
        
        % Calculer le total maximum en respectant le ratio et la tolérance
        min_ratio = target_ratio - tolerance;
        max_ratio = target_ratio + tolerance;
        
        % Essayer d'équilibrer
        for new_total = total_target:-1:1
            test_left = round(new_total * target_ratio);
            test_right = new_total - test_left;
            
            if test_left <= max_left && test_right <= max_right
                test_ratio = test_left / new_total;
                if test_ratio >= min_ratio && test_ratio <= max_ratio
                    n_left = test_left;
                    n_right = test_right;
                    return;
                end
            end
        end
    end
end

function selected_idx = select_cycles_random(cycles_array, n_cycles, seed)
% Sélectionne n_cycles aléatoirement avec un seed spécifique
%
% Entrées :
%   cycles_array : structure array des cycles
%   n_cycles     : nombre de cycles à sélectionner
%   seed         : seed pour le générateur aléatoire
%
% Sortie :
%   selected_idx : indices des cycles sélectionnés (triés)

    total_cycles = length(cycles_array);
    
    if n_cycles >= total_cycles
        selected_idx = 1:total_cycles;
        return;
    end
    
    % Sauvegarder l'état actuel du RNG
    previous_rng = rng;
    
    % Utiliser le seed spécifique
    rng(seed);
    
    % Sélection aléatoire
    selected_idx = randperm(total_cycles, n_cycles);
    selected_idx = sort(selected_idx); % Garder l'ordre chronologique
    
    % Restaurer l'état précédent du RNG
    rng(previous_rng);
end

function summary_table = create_summary_table_balanced(cycle_stats, participants, conditions)
% Crée un tableau récapitulatif incluant les ratios Left/Right
    
    n_participants = length(participants);
    
    % Initialiser les colonnes
    Participant = cell(n_participants, 1);
    Total_Original_Min = zeros(n_participants, 1);
    Total_Selected = zeros(n_participants, 1);
    Selected_Left = zeros(n_participants, 1);
    Selected_Right = zeros(n_participants, 1);
    Left_Ratio = zeros(n_participants, 1);
    
    % Colonnes détaillées par condition
    for iC = 1:length(conditions)
        cond = conditions{iC};
        eval(sprintf('Original_%s_Left = zeros(n_participants, 1);', cond));
        eval(sprintf('Original_%s_Right = zeros(n_participants, 1);', cond));
        eval(sprintf('Original_%s_Total = zeros(n_participants, 1);', cond));
    end
    
    row_idx = 1;
    
    for p = 1:length(participants)
        participant = participants{p};
        
        if ~isfield(cycle_stats, participant)
            continue;
        end
        
        Participant{row_idx} = participant;
        
        % Cycles sélectionnés (commun à toutes conditions)
        Selected_Left(row_idx) = cycle_stats.(participant).selected.Left;
        Selected_Right(row_idx) = cycle_stats.(participant).selected.Right;
        Total_Selected(row_idx) = cycle_stats.(participant).selected.Total;
        Left_Ratio(row_idx) = cycle_stats.(participant).selected.LeftRatio;
        
        % Cycles originaux par condition
        totals = [];
        for iC = 1:length(conditions)
            cond = conditions{iC};
            
            if isfield(cycle_stats.(participant).original, cond)
                left = cycle_stats.(participant).original.(cond).Left;
                right = cycle_stats.(participant).original.(cond).Right;
                total = cycle_stats.(participant).original.(cond).Total;
                
                eval(sprintf('Original_%s_Left(row_idx) = left;', cond));
                eval(sprintf('Original_%s_Right(row_idx) = right;', cond));
                eval(sprintf('Original_%s_Total(row_idx) = total;', cond));
                
                totals = [totals; total];
            end
        end
        
        Total_Original_Min(row_idx) = min(totals);
        
        row_idx = row_idx + 1;
    end
    
    % Créer la table de base
    summary_table = table(Participant(1:row_idx-1), ...
                         Total_Original_Min(1:row_idx-1), ...
                         Total_Selected(1:row_idx-1), ...
                         Selected_Left(1:row_idx-1), ...
                         Selected_Right(1:row_idx-1), ...
                         Left_Ratio(1:row_idx-1), ...
                         'VariableNames', ...
                         {'Participant', 'Total_Original_Min', 'Total_Selected', ...
                          'Selected_Left', 'Selected_Right', 'Left_Ratio_Percent'});
    
    % Convertir le ratio en pourcentage
    summary_table.Left_Ratio_Percent = summary_table.Left_Ratio_Percent * 100;
    
    % Ajouter les colonnes détaillées par condition
    for iC = 1:length(conditions)
        cond = conditions{iC};
        eval(sprintf('summary_table.Original_%s_Left = Original_%s_Left(1:row_idx-1);', cond, cond));
        eval(sprintf('summary_table.Original_%s_Right = Original_%s_Right(1:row_idx-1);', cond, cond));
        eval(sprintf('summary_table.Original_%s_Total = Original_%s_Total(1:row_idx-1);', cond, cond));
    end
end