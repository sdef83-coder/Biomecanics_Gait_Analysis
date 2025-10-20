function fig = radarGaitPlot_5Domains_Inter(SpatioTemporalDATA, condition, groups)
% radarGaitPlot_5Domains_Inter - Radar plot INTER-GROUPES avec 5 domaines
% Compare les 4 groupes d'âge pour une condition donnée

% === Configuration des 5 domaines ===
domains = struct();

% 1. PACE
domains.Pace = {
    'vitFoulee', 'Gait speed', 'Mean';
    'distFoulee', 'Stride length', 'Mean';
    'tempsFoulee', 'Stride time', 'Mean'
};

% 2. RHYTHM
domains.Rhythm = {
    'vitCadencePasParMinute', 'Cadence', 'Mean';
    'NormCadence', 'Normalized cadence', 'Mean';
    'NormWalkRatio', 'Normalized walk ratio', 'Mean'
};

% 3. STABILITY
domains.Stability = {
    'LargeurPas', 'Stride width', 'Mean';
    'MoS_HS_ML', 'MoS HS ML', 'Mean';
    'DoubleSupport', 'Double support', 'Mean'
};

% 4. ASYMMETRY
domains.Asymmetry = {
    'tempsFoulee', 'Stride time', 'SI';
    'distFoulee', 'Stride length', 'SI';
    'LargeurPas', 'Stride width', 'SI'
};

% 5. VARIABILITY
domains.Variability = {
    'tempsFoulee', 'Stride time', 'CV';
    'distFoulee', 'Stride length', 'CV';
    'LargeurPas', 'Stride width', 'CV'
};

% === Collecte des données ===
domainNames = fieldnames(domains);
nDomains = length(domainNames);
nGroups = length(groups);

% Calculer nombre total de variables
nVarsTotal = 0;
for d = 1:nDomains
    nVarsTotal = nVarsTotal + size(domains.(domainNames{d}), 1);
end

% Initialisation
group_means = zeros(nGroups, nVarsTotal);
group_data = cell(nGroups, nVarsTotal);
labels = cell(1, nVarsTotal);
domain_indices = cell(1, nDomains);
stat_type = cell(1, nVarsTotal); % 'mean' ou 'median'

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
        
        % Construire le nom de la variable
        varFullName = [varType '_' varTech];
        
        % Label
        labels{varIdx} = sprintf('%s (%s)', varDisplay, varType);
        
        % Collecter les données pour chaque groupe
        for g = 1:nGroups
            groupName = groups{g};
            
            if isfield(SpatioTemporalDATA, groupName) && ...
               isfield(SpatioTemporalDATA.(groupName), condition)
                
                tableData = SpatioTemporalDATA.(groupName).(condition);
                
                if any(strcmp(tableData.Properties.VariableNames, varFullName))
                    values = tableData.(varFullName);
                    group_data{g, varIdx} = values;
                    
                    % Test de normalité et choix moyenne/médiane
                    if length(values) >= 3 && swtest(values, 0.05) == 1
                        group_means(g, varIdx) = median(values, 'omitnan');
                        stat_type{varIdx} = 'md';
                    else
                        group_means(g, varIdx) = mean(values, 'omitnan');
                        stat_type{varIdx} = 'μ';
                    end
                end
            end
        end
        
        varIdx = varIdx + 1;
    end
end

% === Normalisation ===
all_values = [];
for g = 1:nGroups
    for v = 1:nVarsTotal
        if ~isempty(group_data{g,v})
            all_values = [all_values; group_data{g,v}(:)]; %#ok<AGROW>
        end
    end
end

minVals = zeros(1, nVarsTotal);
maxVals = zeros(1, nVarsTotal);
for v = 1:nVarsTotal
    col_vals = [];
    for g = 1:nGroups
        if ~isempty(group_data{g,v})
            col_vals = [col_vals; group_data{g,v}(:)]; %#ok<AGROW>
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
fig = figure('Name', ['5 Domains Inter - ' condition]);
hold on; axis equal;
set(gca, 'XColor', 'none', 'YColor', 'none');

% === Zones colorées par domaine ===
domainColors = {
    [0.7 0.85 1.0];   % Pace - Bleu clair
    [1.0 0.7 0.7];    % Rhythm - Rouge clair
    [0.7 1.0 0.7];    % Stability - Vert clair
    [1.0 0.85 0.6];   % Asymmetry - Orange clair
    [0.9 0.7 1.0]     % Variability - Violet clair
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
    
    % Ajouter le type de stat au label
    label_with_stat = sprintf('%s (%s)', labels{i}, stat_type{i});
    text(x, y, label_with_stat, 'FontSize', 7, ...
        'HorizontalAlignment', align, 'VerticalAlignment', 'middle');
end

% === Couleurs et marqueurs des groupes ===
groupColors = {
    [0 0.6 1];      % Jeunes Enfants - Bleu
    [1 0.5 0];      % Enfants - Orange
    [0 0.70 0.30];  % Adolescents - Vert
    [0.5 0 1]       % Adultes - Violet
};
groupMarkers = {'^', 'o', 'd', 's'};

% === Tracer les données individuelles (transparentes) ===
for g = 1:nGroups
    color = groupColors{g};
    
    for i = 1:size(group_data{g,1}, 1)
        r_individual = zeros(1, nVarsTotal);
        for v = 1:nVarsTotal
            if ~isempty(group_data{g,v}) && i <= length(group_data{g,v})
                r_individual(v) = scale(group_data{g,v}(i), v);
            end
        end
        r_individual = [r_individual, r_individual(1)];
        
        plot(r_individual .* cos(theta), r_individual .* sin(theta), ...
            '-', 'Color', [color 0.2], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end

% === Tracer les moyennes/médianes (lignes épaisses) ===
for g = 1:nGroups
    groupName = groups{g};
    color = groupColors{g};
    marker = groupMarkers{g};
    
    % Normaliser les moyennes
    r_norm = arrayfun(@(j) scale(group_means(g,j), j), 1:nVarsTotal);
    r_norm = [r_norm, r_norm(1)];
    
    x = r_norm .* cos(theta);
    y = r_norm .* sin(theta);
    
    plot(x, y, '-', 'Color', color, 'LineWidth', 2.5, ...
        'Marker', marker, 'MarkerSize', 8, 'MarkerFaceColor', color, ...
        'DisplayName', groupName);
    
    % Valeurs numériques
    for i = 1:nVarsTotal
        text(x(i)*1.10, y(i)*1.10, sprintf('%.2f', group_means(g,i)), ...
            'FontSize', 6, 'Color', color, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold');
    end
end

% === Légendes ===
legend('Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 10);

% Légende des domaines
annotation('textbox', [0.01 0.60 0.22 0.30], 'String', {
    '\bfDomains:', ...
    '\color[rgb]{0.35 0.6 1.0}█ Pace', ...
    '\color[rgb]{1.0 0.3 0.3}█ Rhythm', ...
    '\color[rgb]{0.3 0.8 0.3}█ Stability', ...
    '\color[rgb]{1.0 0.6 0.3}█ Asymmetry', ...
    '\color[rgb]{0.7 0.5 1.0}█ Variability'}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'EdgeColor', 'black', 'FontSize', 9, 'LineWidth', 1);

% Légende statistiques
annotation('textbox', [0.01 0.40 0.22 0.15], 'String', {
    '\bfStatistics:', ...
    '(μ) = Mean used', ...
    '(md) = Median used'}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'EdgeColor', 'black', 'FontSize', 9, 'LineWidth', 1);

% Titre
annotation('textbox', [0 0.955 1 0.05], ...
    'String', ['\bf5 Gait Domains - ' condition ' - Inter-groups comparison'], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FontSize', 14, 'FontWeight', 'bold');

end