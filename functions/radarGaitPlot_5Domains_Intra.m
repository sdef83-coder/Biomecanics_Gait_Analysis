function fig = radarGaitPlot_5Domains_Intra(SpatioTemporalDATA, group, conditions, condColors)
% radarGaitPlot_5Domains_Intra - Radar plot INTRA-GROUPE avec 5 domaines
% Compare les 3 conditions pour un groupe d'âge donné

% === Configuration des 5 domaines ===
domains = struct();

domains.Pace = {
    'vitFoulee', 'Gait speed', 'Mean';
    'distFoulee', 'Stride length', 'Mean';
    'tempsFoulee', 'Stride time', 'Mean'
};

domains.Rhythm = {
    'vitCadencePasParMinute', 'Cadence', 'Mean';
    'NormCadence', 'Normalized cadence', 'Mean';
    'NormWalkRatio', 'Normalized walk ratio', 'Mean'
};

domains.Stability = {
    'LargeurPas', 'Stride width', 'Mean';
    'MoS_HS_ML', 'MoS HS ML', 'Mean';
    'DoubleSupport', 'Double support', 'Mean'
};

domains.Asymmetry = {
    'tempsFoulee', 'Stride time', 'SI';
    'distFoulee', 'Stride length', 'SI';
    'LargeurPas', 'Stride width', 'SI'
};

domains.Variability = {
    'tempsFoulee', 'Stride time', 'CV';
    'distFoulee', 'Stride length', 'CV';
    'LargeurPas', 'Stride width', 'CV'
};

% === Collecte des données ===
domainNames = fieldnames(domains);
nDomains = length(domainNames);
nConds = length(conditions);

nVarsTotal = 0;
for d = 1:nDomains
    nVarsTotal = nVarsTotal + size(domains.(domainNames{d}), 1);
end

% Initialisation
cond_means = zeros(nConds, nVarsTotal);
cond_data = cell(nConds, nVarsTotal);
labels = cell(1, nVarsTotal);
domain_indices = cell(1, nDomains);
stat_type = cell(1, nVarsTotal);

varIdx = 1;
for d = 1:nDomains
    domainName = domainNames{d};
    domainVars = domains.(domainName);
    nVars = size(domainVars, 1);
    
    domain_indices{d} = varIdx:(varIdx + nVars - 1);
    
    for v = 1:nVars
        varTech = domainVars{v, 1};
        varDisplay = domainVars{v, 2};
        varType = domainVars{v, 3};
        
        varFullName = [varType '_' varTech];
        labels{varIdx} = sprintf('%s (%s)', varDisplay, varType);
        
        for c = 1:nConds
            cond = conditions{c};
            
            if isfield(SpatioTemporalDATA, group) && ...
               isfield(SpatioTemporalDATA.(group), cond)
                
                tableData = SpatioTemporalDATA.(group).(cond);
                
                if any(strcmp(tableData.Properties.VariableNames, varFullName))
                    values = tableData.(varFullName);
                    cond_data{c, varIdx} = values;
                    
                    if length(values) >= 3 && swtest(values, 0.05) == 1
                        cond_means(c, varIdx) = median(values, 'omitnan');
                        stat_type{varIdx} = 'md';
                    else
                        cond_means(c, varIdx) = mean(values, 'omitnan');
                        stat_type{varIdx} = 'μ';
                    end
                end
            end
        end
        
        varIdx = varIdx + 1;
    end
end

% === Normalisation ===
minVals = zeros(1, nVarsTotal);
maxVals = zeros(1, nVarsTotal);
for v = 1:nVarsTotal
    col_vals = [];
    for c = 1:nConds
        if ~isempty(cond_data{c,v})
            col_vals = [col_vals; cond_data{c,v}(:)]; %#ok<AGROW>
        end
    end
    if ~isempty(col_vals)
        minVals(v) = min(col_vals);
        maxVals(v) = max(col_vals);
    end
end

scale = @(val, j) (val - minVals(j)) / max(maxVals(j) - minVals(j), eps);

% === Angles ===
theta = linspace(0, 2*pi, nVarsTotal+1);
theta(end) = [];
theta = [theta, theta(1)];

% === Figure ===
fig = figure('Name', ['5 Domains Intra - ' group]);
hold on; axis equal;
set(gca, 'XColor', 'none', 'YColor', 'none');

% === Zones colorées ===
domainColors = {
    [0.7 0.85 1.0];
    [1.0 0.7 0.7];
    [0.7 1.0 0.7];
    [1.0 0.85 0.6];
    [0.9 0.7 1.0]
};

for d = 1:nDomains
    idx = domain_indices{d};
    col = domainColors{d};
    
    for i = idx
        if i < length(theta)
            patch_theta = [theta(i), theta(i+1), 0];
            patch_r = [1 1 0];
            [x, y] = pol2cart(patch_theta, patch_r);
            fill(x, y, col, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
end

% === Axes et labels ===
for i = 1:nVarsTotal
    plot([0 cos(theta(i))], [0 sin(theta(i))], 'k--', 'HandleVisibility', 'off');
    
    ang = rad2deg(theta(i));
    x = 1.18 * cos(theta(i));
    y = 1.18 * sin(theta(i));
    
    align = 'left';
    if ang > 90 && ang < 270
        align = 'right';
    end
    
    label_with_stat = sprintf('%s (%s)', labels{i}, stat_type{i});
    text(x, y, label_with_stat, 'FontSize', 7, ...
        'HorizontalAlignment', align, 'VerticalAlignment', 'middle');
end

% === Tracer données individuelles ===
for c = 1:nConds
    color = condColors{c};
    
    for i = 1:size(cond_data{c,1}, 1)
        r_individual = zeros(1, nVarsTotal);
        for v = 1:nVarsTotal
            if ~isempty(cond_data{c,v}) && i <= length(cond_data{c,v})
                r_individual(v) = scale(cond_data{c,v}(i), v);
            end
        end
        r_individual = [r_individual, r_individual(1)];
        
        plot(r_individual .* cos(theta), r_individual .* sin(theta), ...
            '-', 'Color', [color 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end

% === Tracer moyennes ===
for c = 1:nConds
    cond = conditions{c};
    color = condColors{c};
    
    r_norm = arrayfun(@(j) scale(cond_means(c,j), j), 1:nVarsTotal);
    r_norm = [r_norm, r_norm(1)];
    
    x = r_norm .* cos(theta);
    y = r_norm .* sin(theta);
    
    plot(x, y, '-', 'Color', color, 'LineWidth', 2.5, ...
        'Marker', 'o', 'MarkerSize', 8, 'MarkerFaceColor', color, ...
        'DisplayName', cond);
    
    for i = 1:nVarsTotal
        text(x(i)*1.10, y(i)*1.10, sprintf('%.2f', cond_means(c,i)), ...
            'FontSize', 6, 'Color', color, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold');
    end
end

% === Légendes ===
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 10);

annotation('textbox', [0.01 0.60 0.22 0.30], 'String', {
    '\bfDomains:', ...
    '\color[rgb]{0.35 0.6 1.0}█ Pace', ...
    '\color[rgb]{1.0 0.3 0.3}█ Rhythm', ...
    '\color[rgb]{0.3 0.8 0.3}█ Stability', ...
    '\color[rgb]{1.0 0.6 0.3}█ Asymmetry', ...
    '\color[rgb]{0.7 0.5 1.0}█ Variability'}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'EdgeColor', 'black', 'FontSize', 9, 'LineWidth', 1);

annotation('textbox', [0.01 0.40 0.22 0.15], 'String', {
    '\bfStatistics:', ...
    '(μ) = Mean', ...
    '(md) = Median'}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'EdgeColor', 'black', 'FontSize', 9, 'LineWidth', 1);

annotation('textbox', [0 0.955 1 0.05], ...
    'String', ['\bf5 Gait Domains - ' group ' - Intra-group comparison'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FontSize', 14, 'FontWeight', 'bold');

end