function [results] = nonParametricANOVA(data, groupVar, conditionVar, participantVar)
    results = struct();
    alpha = 0.05;

    % Test Kruskal-Wallis : effet du groupe
    try
        [p_group, ~, stats_group] = kruskalwallis(data, groupVar, 'off');
        results.group_effect.test = 'Kruskal-Wallis';
        results.group_effect.p_value = p_group;
        results.group_effect.significant = p_group < alpha;
        results.group_effect.stats = stats_group;

        if p_group < alpha
            c_group = multcompare(stats_group, 'Display', 'off');
            results.group_effect.posthoc = c_group;
        end
    catch
        warning('❌ Erreur Kruskal-Wallis');
        results.group_effect.test = 'Failed';
        results.group_effect.p_value = NaN;
        results.group_effect.significant = false;
    end

    % Test Friedman : effet des conditions (mesures répétées)
    try
        uniqueSubjects = unique(participantVar);
        uniqueConditions = unique(conditionVar);

        friedmanData = NaN(length(uniqueSubjects), length(uniqueConditions));

        for i = 1:length(uniqueSubjects)
            subjIdx = participantVar == uniqueSubjects(i);
            subjData = data(subjIdx);
            subjCond = conditionVar(subjIdx);

            for j = 1:length(uniqueConditions)
                condIdx = subjCond == uniqueConditions(j);
                if any(condIdx)
                    friedmanData(i, j) = subjData(condIdx);
                end
            end
        end

        [p_friedman, ~, stats_friedman] = friedman(friedmanData, 1, 'off');
        results.condition_effect.test = 'Friedman';
        results.condition_effect.p_value = p_friedman;
        results.condition_effect.significant = p_friedman < alpha;
        results.condition_effect.stats = stats_friedman;

    catch
        warning('❌ Erreur Friedman');
        results.condition_effect.test = 'Failed';
        results.condition_effect.p_value = NaN;
        results.condition_effect.significant = false;
    end
end