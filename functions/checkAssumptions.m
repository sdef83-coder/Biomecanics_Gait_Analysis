function [isNormal, isHomogeneous, results] = checkAssumptions(data, groupVar, conditionVar)
    results = struct();
    alpha = 0.05;

    % Supprimer les NaN
    validIdx = ~isnan(data) & ~isnan(double(groupVar)) & ~isnan(double(conditionVar));
    data = data(validIdx);
    groupVar = groupVar(validIdx);
    conditionVar = conditionVar(validIdx);

    % 1. Normalité des résidus via Shapiro-Wilk
    try
        [~, ~, stats] = anova1(data(:), groupVar, 'off');
        residuals = stats.resid;

        [h_sw, p_sw] = swtest(residuals);
results.normality.test = 'Shapiro-Wilk (swtest)';
results.normality.p_value = p_sw;
results.normality.is_normal = ~h_sw;

        isNormal = results.normality.is_normal;
    catch
        warning('❌ Erreur test de normalité (Shapiro-Wilk)');
        isNormal = false;
        results.normality.test = 'Failed';
        results.normality.p_value = NaN;
        results.normality.is_normal = false;
    end

    % 2. Homogénéité des variances (Levene)
    try
        [p_levene, ~] = vartestn(data, groupVar, 'TestType', 'LeveneAbsolute', 'Display', 'off');
        results.homogeneity.test = 'Levene';
        results.homogeneity.p_value = p_levene;
        results.homogeneity.is_homogeneous = p_levene > alpha;

        isHomogeneous = results.homogeneity.is_homogeneous;
    catch
        warning('❌ Erreur test Levene');
        isHomogeneous = false;
        results.homogeneity.test = 'Failed';
        results.homogeneity.p_value = NaN;
        results.homogeneity.is_homogeneous = false;
    end

    % 3. Sphéricité (approximation via normalité des différences)
    try
        uniqueConditions = unique(conditionVar);
        results.sphericity.test = 'Sphericité (différences)';
        results.sphericity.p_value = NaN;
        results.sphericity.is_spherical = true;

        if length(uniqueConditions) > 2
            differences = [];
            for i = 1:length(uniqueConditions)-1
                for j = i+1:length(uniqueConditions)
                    cond1_idx = conditionVar == uniqueConditions(i);
                    cond2_idx = conditionVar == uniqueConditions(j);

                    if sum(cond1_idx) == sum(cond2_idx)
                        diff_values = data(cond1_idx) - data(cond2_idx);
                        differences = [differences; diff_values];
                    end
                end
            end

            if ~isempty(differences)
                [h_sph, p_sph] = swtest(differences);
                results.sphericity.p_value = p_sph;
                results.sphericity.is_spherical = ~h_sph;
            end
        end
    end