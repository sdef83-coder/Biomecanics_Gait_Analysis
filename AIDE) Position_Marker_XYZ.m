%% export_STRN_filtered_trials.m
% Exporte les positions du marqueur sternum (STRN) filtrées (LP 6 Hz, Butter ordre 2, filtfilt)
% pour qu'un tiers puisse recalculer SPARC / LDLJ / jerk, etc... indépendamment et valider la méthode.
%
% Dépendances: BTK (btkReadAcquisition, btkGetMarkers, btkCloseAcquisition)
%
% Sorties:
%  - 1 fichier .mat global (toutes surfaces + essais)
%  - 1 fichier .csv par essai (optionnel, recommandé)

clear; clc; close all;

%% ===================== PARAMETRES UTILISATEUR =====================
% BTK
addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\Scripts\btk'));

% Identifiants / structure des fichiers
sujet_id  = 'CTL_33';
surfaces  = {'Plat','Medium','High'};
essais    = 1:10;

% Dossier contenant les C3D (adapter)
base_dir  = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\adults';

% Acquisition / filtrage (aligné sur ton script)
fs        = 100;   % Hz
fc        = 6;     % Hz
filter_order = 2;

% Export
export_dir = 'C:\Users\silve\Desktop\EXPORT_STRN_for_Supervisor';
export_csv_per_trial = true;  % true = 1 CSV par essai
save_raw_too = true;          % true = on sauvegarde STRN brut + filtré

%% ===================== INITIALISATION =====================
if ~exist(export_dir, 'dir')
    mkdir(export_dir);
end

[b,a] = butter(filter_order, fc/(fs/2), 'low');

S = struct();
S.meta.sujet_id = sujet_id;
S.meta.surfaces = surfaces;
S.meta.essais   = essais;
S.meta.fs       = fs;
S.meta.fc       = fc;
S.meta.filter_order = filter_order;
S.meta.marker_name  = 'STRN';
S.data = [];   % sera un struct array par essai

log_errors = {};
k = 0;

fprintf('=== Export STRN filtré | Sujet %s ===\n', sujet_id);

%% ===================== BOUCLE SUR FICHIERS =====================
for s = 1:numel(surfaces)
    surface = surfaces{s};

    for essai = essais
        filename = sprintf('%s_%s_%02d.c3d', sujet_id, surface, essai);
        c3d_path = fullfile(base_dir, filename);

        if ~isfile(c3d_path)
            msg = sprintf('Fichier manquant: %s', filename);
            log_errors{end+1} = msg;
            fprintf('  [MISS] %s\n', msg);
            continue;
        end

        try
            acq = btkReadAcquisition(c3d_path);
            markers = btkGetMarkers(acq);

            if ~isfield(markers, 'STRN') || isempty(markers.STRN) || size(markers.STRN,2) ~= 3
                msg = sprintf('STRN absent/invalide: %s', filename);
                log_errors{end+1} = msg;
                fprintf('  [SKIP] %s\n', msg);
                btkCloseAcquisition(acq);
                continue;
            end

            STRN_raw = markers.STRN; % [N x 3]
            N = size(STRN_raw,1);
            t = (0:N-1)'/fs;

            % Filtrage (identique à ta logique: filtfilt sur positions)
            STRN_filt = filtfilt(b, a, STRN_raw);

            % Stockage
            k = k + 1;
            S.data(k).filename   = filename;
            S.data(k).surface    = surface;
            S.data(k).essai      = essai;
            S.data(k).fs         = fs;
            S.data(k).t          = t;
            S.data(k).STRN_filt  = STRN_filt;

            if save_raw_too
                S.data(k).STRN_raw = STRN_raw;
            end

            % Export CSV par essai
            if export_csv_per_trial
                T = table(t, STRN_filt(:,1), STRN_filt(:,2), STRN_filt(:,3), ...
                    'VariableNames', {'time_s','STRN_X_filt','STRN_Y_filt','STRN_Z_filt'});

                if save_raw_too
                    T.STRN_X_raw = STRN_raw(:,1);
                    T.STRN_Y_raw = STRN_raw(:,2);
                    T.STRN_Z_raw = STRN_raw(:,3);
                end

                csv_name = sprintf('%s_%s_%02d_STRN_filtered.csv', sujet_id, surface, essai);
                writetable(T, fullfile(export_dir, csv_name));
            end

            btkCloseAcquisition(acq);

            fprintf('  [OK] %s | N=%d frames\n', filename, N);

        catch ME
            msg = sprintf('Erreur %s : %s', filename, ME.message);
            log_errors{end+1} = msg;
            fprintf('  [ERR] %s\n', msg);
            try, btkCloseAcquisition(acq); catch, end %#ok<CTCH>
        end
    end
end

%% ===================== SAUVEGARDE MAT GLOBAL =====================
mat_name = sprintf('%s_STRN_filtered_allTrials.mat', sujet_id);
save(fullfile(export_dir, mat_name), 'S', 'log_errors', '-v7.3');

fprintf('\n=== Terminé ===\n');
fprintf('Essais exportés: %d\n', numel(S.data));
fprintf('Dossier export : %s\n', export_dir);

if ~isempty(log_errors)
    fprintf('Avertissements/erreurs: %d (voir log_errors dans .mat)\n', numel(log_errors));
end