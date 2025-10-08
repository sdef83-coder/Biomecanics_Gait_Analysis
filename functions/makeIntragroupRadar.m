function fig = makeIntragroupRadar(DATA, group, cond_list, param_names, colors)

nCond = length(cond_list);
nParams = length(param_names);
nVar = nParams * 2;

rename_vars = containers.Map( ...
    {'pctSimpleAppuie','DoubleSupport','LargeurPas','vitFoulee','distFoulee','NormWalkRatio','vitCadencePasParMinute', 'NormStepLength', 'NormCadence'}, ...
    {'Single support time (%)','Double support time (%)','Stride width (mm)', ...
     'Gait speed (m.s^{-1})','Stride length (m)', ...
     'Normalized Walk ratio (ua)','Cadence (step.min^{-1})', 'Normalized Step length (ua)', 'Normalized Cadence (ua)'});

% === Initialisation
data = cell(1, nCond);
means = zeros(nCond, nVar);
stat_labels = repmat({''}, 1, nVar);  % Stocke μ ou md par variable

% === Extraction + moyenne/médiane
for c = 1:nCond
    cond = cond_list{c};
    temp = zeros(length(DATA.(group).(cond).(['Mean_' param_names{1}])), nVar);
    for i = 1:nParams
        base = param_names{i};
        temp(:,i) = DATA.(group).(cond).(['Mean_' base]);
        temp(:,i+nParams) = DATA.(group).(cond).(['CV_' base]);
    end
    data{c} = temp;

    for j = 1:nVar
        if swtest(temp(:,j), 0.05) == 1
            means(c,j) = median(temp(:,j), 'omitnan');
            stat_labels{j} = '(md)';
        else
            means(c,j) = mean(temp(:,j), 'omitnan');
            % si ce n'est pas déjà "md" par une autre condition
            if ~strcmp(stat_labels{j}, '(md)')
                stat_labels{j} = '(μ)';
            end
        end
    end
end

% === Création des labels avec (μ)/(md)
labels = cell(1, nVar);
for i = 1:nParams
    base = param_names{i};
    if isKey(rename_vars, base)
        nom = rename_vars(base);
        labels{i} = [nom ' ' stat_labels{i}];
        labels{i+nParams} = [extractBefore(nom, ' (') ' (CV) ' stat_labels{i+nParams}];
    else
        labels{i} = [base ' ' stat_labels{i}];
        labels{i+nParams} = [base ' (CV) ' stat_labels{i+nParams}];
    end
end

% === Normalisation
all_data = cat(1, data{:});
minVals = min(all_data,[],1);
maxVals = max(all_data,[],1);
scale = @(val,j) (val - minVals(j)) / (maxVals(j) - minVals(j));

theta = linspace(0, 2*pi, nVar+1); theta(end) = []; theta = [theta, theta(1)];

% === Figure
fig = figure('Name',['Intragroupe - ' group]); hold on; axis equal;
set(gca, 'XColor', 'none', 'YColor', 'none');

% === Zones thématiques
zones = {
    [1 2 3], [0.7 0.85 1.0];
    [4 5 6], [1.0 0.7 0.7];
    [7 8 9], [0.7 1.0 0.7];
    10:18, [1.0 1.0 0.6];
};
for z = 1:size(zones,1)
    idx = zones{z,1}; col = zones{z,2};
    for i = idx
        patch_theta = [theta(i), theta(i+1), 0];
        patch_r = [1 1 0];
        [x, y] = pol2cart(patch_theta, patch_r);
        fill(x, y, col, 'FaceAlpha', 0.25, 'EdgeColor','none', 'HandleVisibility','off');
    end
end

% === Tracer axes + labels
for i = 1:nVar
    plot([0 cos(theta(i))], [0 sin(theta(i))], 'k--', 'HandleVisibility', 'Off');
    ang = rad2deg(theta(i));
    x = 1.05 * cos(theta(i)); y = 1.05 * sin(theta(i));
    align = 'left'; if ang > 90 && ang < 270, align = 'right'; end
    text(x, y, labels{i}, 'FontSize', 8, 'HorizontalAlignment', align, 'VerticalAlignment','middle');
end

% === Tracer les courbes
for c = 1:nCond
    r_mean = arrayfun(@(j) scale(means(c,j), j), 1:nVar);
    r_mean = [r_mean r_mean(1)];
    x = r_mean .* cos(theta); y = r_mean .* sin(theta);
    plot(x, y, '-', 'Color', colors{c}, 'LineWidth', 2, 'DisplayName', cond_list{c});

    for i = 1:nVar
        text(x(i)*1.05, y(i)*1.05, sprintf('%.2f', means(c,i)), ...
            'FontSize', 8, 'Color', 'k', 'HorizontalAlignment','center');
    end

    for i = 1:size(data{c},1)
        r = arrayfun(@(j) scale(data{c}(i,j), j), 1:nVar);
        r = [r r(1)];
        plot(r .* cos(theta), r .* sin(theta), '-', 'Color', [colors{c} 0.3], 'HandleVisibility', 'Off');
    end
end

% === Légende des conditions
legend('Location','southoutside','Orientation','horizontal');

% === Légende des statistiques (au-dessus des zones thématiques)
annotation('textbox', [0.09 0.32 0.2 0.05], 'String', {
    '\bfStatistiques :', ...
    '(μ) : Moyenne utilisée', ...
    '(md) : Médiane utilisée'}, ...
    'FitBoxToText','on', ...
    'EdgeColor','none', ...
    'FontSize', 9);

% === Légende zones thématiques
annotation('textbox', [0.09 0.10 0.2 0.15], 'String', ...
    {'\color[rgb]{0.35 0.6 1.0} Pace', ...
     '\color[rgb]{1.0 0.3 0.3} Rhythm', ...
     '\color[rgb]{0.3 0.8 0.3} Stability', ...
     '\color[rgb]{0.85 0.85 0} Variability'}, ...
    'FitBoxToText','on', 'BackgroundColor','none', ...
    'EdgeColor','none', 'FontSize', 9);

% === Titre
annotation('textbox', [0 0.955 1 0.05], ...
    'String', ['\bfSpatiotemporal variables - ' group], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FontSize', 14);

end