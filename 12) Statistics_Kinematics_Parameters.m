%% ENFANTS VS ADULTES
clc; clear; close all;

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\gaitAnalysisGUI\result\Matrice');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\PROGRAMMATION\spm1dmatlab-master'));

Enfants = load("MATRICES_Enfants.mat");
Adultes = load('MATRICES_Adultes.mat');
JeunesEnfants = load("MATRICES_JeunesEnfants.mat");

groupList = {
    struct('G1', Enfants, 'G1name', 'Enfants', 'color1', 'r', ...
           'G2', Adultes, 'G2name', 'Adultes', 'color2', 'b'), ...
    struct('G1', Enfants, 'G1name', 'Enfants', 'color1', 'r', ...
           'G2', JeunesEnfants, 'G2name', 'JeunesEnfants', 'color2', 'm'), ...
    struct('G1', JeunesEnfants, 'G1name', 'JeunesEnfants', 'color1', 'm', ...
           'G2', Adultes, 'G2name', 'Adultes', 'color2', 'b')
};

saveFolder = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\gaitAnalysisGUI\result\Fig';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder)
end

x = 1:100;

% Fonction affichage bande ±1 SD
affiche_avec_std = @(mean_data, std_data, color) fill([x fliplr(x)], ...
    [mean_data + std_data, fliplr(mean_data - std_data)], ...
    color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Listes
plans = {'Sagittal', 'Frontal', 'Transverse'};
nomPlans = {'Plan Sagittal', 'Plan Frontal', 'Plan Transverse'};
articulations = {'Ankle', 'Knee', 'Hip'};
titres = {'Cheville', 'Genou', 'Hanche'};
conditions = {'Plat', 'Medium', 'High'};
nConditions = numel(conditions);

for g = 1:length(groupList)
    grp = groupList{g};

    for p = 1:length(plans)
        plan = plans{p};
        nomPlan = nomPlans{p};

        for a = 1:length(articulations)
            art = articulations{a};
            nomArt = titres{a};

            fig = figure('Position',[100 100 1400 800]);
            sgtitle(['Comparaison ' grp.G1name ' vs ' grp.G2name ' - ' nomPlan ' - ' nomArt]);

            for c = 1:nConditions
                cond = conditions{c};

                subplot(nConditions,2,(c-1)*2+1)
                hold on
                mean_1 = grp.G1.MATRICES.(plan).(art).(['Mean_' cond]);
                std_1  = grp.G1.MATRICES.(plan).(art).(['Std_' cond]);
                mean_2 = grp.G2.MATRICES.(plan).(art).(['Mean_' cond]);
                std_2  = grp.G2.MATRICES.(plan).(art).(['Std_' cond]);

                affiche_avec_std(mean_1, std_1, grp.color1)
                affiche_avec_std(mean_2, std_2, grp.color2)

                plot(x, mean_1, grp.color1, 'LineWidth', 1.5)
                plot(x, mean_2, grp.color2, 'LineWidth', 1.5)

                title(['Cinématique - ' cond])
                legend(grp.G1name, grp.G2name)
                xlabel('% cycle de marche')
                ylabel('Angle (°)')
                grid on

                subplot(nConditions,2,(c-1)*2+2)
                hold on
                Y1 = grp.G1.MATRICES.(plan).(art).([cond]);
                Y2 = grp.G2.MATRICES.(plan).(art).([cond]);

                valid_idx_Y1 = all(~isnan(Y1), 2);
                valid_idx_Y2 = all(~isnan(Y2), 2);
                Y1 = Y1(valid_idx_Y1, :);
                Y2 = Y2(valid_idx_Y2, :);

                alpha = 0.05;
                normtest = spm1d.stats.normality.ttest2(Y1, Y2);
                normi = normtest.inference(alpha, 'interp', true);
                non_normal_points = sum(normi.h0reject);

                if non_normal_points > 0.05 * length(normi.h0reject)
                    disp(['⚠️  SnPM utilisé pour ' plan ' - ' art ' - ' cond ...
                        ' : normalité violée à ' num2str(non_normal_points) ' points']);
                    spm = spm1d.stats.nonparam.ttest2(Y1, Y2);
                    nPermMax = spm.permuter.nPermTotal;  % nombre max possible de permutations
                    iterations = min(5000, nPermMax);   % 5000 si possible, sinon le max
                    spmi = spm.inference(alpha, 'two_tailed', true, 'interp', true, 'iterations', iterations);

                else
                    disp(['✅ SPM utilisé pour ' plan ' - ' art ' - ' cond]);
                    spm = spm1d.stats.ttest2(Y1, Y2);
                    spmi = spm.inference(alpha, 'two_tailed', true, 'interp', true);
                end

                spmi.plot();
                spmi.plot_threshold_label();
                spmi.plot_p_values(); 
                title(['SPM1D - ' cond])
                ylabel('T-statistique')
                xlabel('% cycle de marche')
                grid on

                if spmi.h0reject && isstruct(spmi.clusters)
                    for i = 1:length(spmi.clusters)
                        if isfield(spmi.clusters(i), 'h0reject') && spmi.clusters(i).h0reject
                            cluster = spmi.clusters(i).indices;
                            fill([x(cluster) fliplr(x(cluster))], ...
                                 [max(ylim)*ones(size(cluster)) fliplr(min(ylim)*ones(size(cluster)))], ...
                                 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                            text(mean(x(cluster)), max(ylim)*0.9, '*', ...
                                'FontSize', 16, 'Color', 'k', 'FontWeight', 'bold', ...
                                'HorizontalAlignment', 'center');
                        end
                    end
                else
                    disp(['→ Aucune différence significative pour ' plan ' - ' art ' - ' cond]);
                end
            end

            % Sauvegarde de la figure
            saveName = fullfile(saveFolder, ['Comparaison_' grp.G1name '_vs_' grp.G2name '_' nomArt '_' plan]);
            exportgraphics(fig, [saveName '.png'], 'Resolution', 300)
            close(fig)
        end
    end
end
%%
groupes = {'Enfants', 'Adultes', 'JeunesEnfants'};
colors = struct();
colors.Enfants = {'r', 'm'};
colors.Adultes = {'b', 'c'};
colors.JeunesEnfants = {'g', [0.5 0.8 0.2]};

for p = 1:length(plans)
    plan = plans{p};
    nomPlan = nomPlans{p};

    for a = 1:length(articulations)
        art = articulations{a};
        nomArt = titres{a};

        for g = 1:length(groupes)
            grp = groupes{g};
            DATA = eval(grp);  % charge 'Enfants' ou 'Adultes'
            fig = figure('Position', [100 100 1400 400]);
            sgtitle(['Comparaison intra-groupe Plat vs High - ' grp ' - ' nomPlan ' - ' nomArt]);

            % ==== GRAPHIQUE CINÉMATIQUE
            subplot(1,2,1); hold on;
            mean_plat = DATA.MATRICES.(plan).(art).Mean_Plat;
            std_plat  = DATA.MATRICES.(plan).(art).Std_Plat;
            mean_high = DATA.MATRICES.(plan).(art).Mean_High;
            std_high  = DATA.MATRICES.(plan).(art).Std_High;

            affiche_avec_std(mean_plat, std_plat, colors.(grp){1});
            affiche_avec_std(mean_high, std_high, colors.(grp){2});
            plot(x, mean_plat, colors.(grp){1}, 'LineWidth', 1.5);
            plot(x, mean_high, colors.(grp){2}, 'LineWidth', 1.5);

            legend('Plat', 'High');
            title('Cinématique Plat vs High');
            xlabel('% cycle de marche');
            ylabel('Angle (°)');
            grid on;

            % ==== ANALYSE SPM
            subplot(1,2,2); hold on;
            Y1 = DATA.MATRICES.(plan).(art).Plat;
            Y2 = DATA.MATRICES.(plan).(art).High;
            Y1 = Y1(all(~isnan(Y1), 2), :);
            Y2 = Y2(all(~isnan(Y2), 2), :);

            alpha = 0.05;
            normtest = spm1d.stats.normality.ttest2(Y1, Y2);
            normi = normtest.inference(alpha, 'interp', true);
            non_normal_points = sum(normi.h0reject);

            if non_normal_points > 0.05 * length(normi.h0reject)
                disp(['⚠️  SnPM utilisé pour intra-groupe ' grp ' - ' plan ' - ' art]);
                spm = spm1d.stats.nonparam.ttest2(Y1, Y2);
                spmi = spm.inference(alpha, 'two_tailed', true, 'interp', true, 'iterations', 5000);
            else
                disp(['✅ SPM utilisé pour intra-groupe ' grp ' - ' plan ' - ' art]);
                spm = spm1d.stats.ttest2(Y1, Y2);
                spmi = spm.inference(alpha, 'two_tailed', true, 'interp', true);
            end

            spmi.plot();
            spmi.plot_threshold_label();
            spmi.plot_p_values(); 
            title('SPM1D Plat vs High');
            ylabel('T-statistique');
            xlabel('% cycle de marche');
            grid on;

            if spmi.h0reject && isstruct(spmi.clusters)
                for i = 1:length(spmi.clusters)
                    if isfield(spmi.clusters(i), 'h0reject') && spmi.clusters(i).h0reject
                        cluster = spmi.clusters(i).indices;
                        fill([x(cluster) fliplr(x(cluster))], ...
                            [max(ylim)*ones(size(cluster)) fliplr(min(ylim)*ones(size(cluster)))], ...
                            'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                        text(mean(x(cluster)), max(ylim)*0.9, '*', ...
                            'FontSize', 16, 'Color', 'k', 'FontWeight', 'bold', ...
                            'HorizontalAlignment', 'center');
                    end
                end
            else
                disp(['→ Aucune différence significative intra-groupe pour ' grp ' - ' plan ' - ' art]);
            end

            % ==== SAUVEGARDE
            saveName = fullfile(saveFolder, ['Intragroupe_' grp '_' nomArt '_' plan]);
            exportgraphics(fig, [saveName '.png'], 'Resolution', 300);
            close(fig);
        end
    end
end