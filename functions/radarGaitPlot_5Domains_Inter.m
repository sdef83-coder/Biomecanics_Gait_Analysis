function fig = radarGaitPlot_5Domains_Inter(SpatioTemporalDATA, condition, groups)
% Radar plot INTER-GROUPES avec 5 domaines

% ====== paramètres d’affichage ======
LABEL_RADIUS_FACT = 1.05;                 % distance des labels
SHOW_NUM_VALUES   = true;                 % afficher valeurs numériques
FIG_POS           = [100 80 1500 950];    % taille de la figure

% === 0. Mapping "technique" -> "lisible" ===
nameMap = containers.Map;
nameMap('vitFoulee')               = 'Gait speed (m.s^{-1})';
nameMap('distFoulee')              = 'Stride length (m)';
nameMap('NormStepLength')          = 'Norm Step length (ua)';
nameMap('tempsFoulee')             = 'Stride time (s)';
nameMap('vitCadencePasParMinute')  = 'Cadence (step.min^{-1})';
nameMap('NormCadence')             = 'Norm Cadence (ua)';
nameMap('NormWalkRatio')           = 'Norm WR (ua)';
nameMap('LargeurPas')              = 'Stride width (cm)';
nameMap('DoubleSupport')           = 'Double support time (%)';
% MoS normalisées (%L0)
nameMap('MoS_AP_HS_pL0')           = 'MoS AP HS (%L0)';
nameMap('MoS_ML_HS_pL0')           = 'MoS ML HS (%L0)';

% === 1. Domaines ===
domains = struct();

domains.Pace = {
    'vitFoulee',      'Gait speed (m/s)',   'Mean';
    'NormStepLength', 'Norm Stride length (ua)', 'Mean';
    'tempsFoulee',    'Stride time(s)',          'Mean'
};

domains.Rhythm = {
    'vitCadencePasParMinute', 'Cadence (step/min)',   'Mean';
    'NormCadence',            'Norm Cadence (ua)',         'Mean';
    'NormWalkRatio',          'Norm Walk ratio (ua)',      'Mean'
};

domains.Stability = {
    'LargeurPas',     'Stride width (cm)',      'Mean';
    'MoS_ML_HS_pL0',  'MoS HS ML (%L0)',        'Mean';
    'MoS_AP_HS_pL0',  'MoS HS AP (%L0)',        'Mean';
    'DoubleSupport',  'Double support (%)',     'Mean'
};

domains.Asymmetry = {
    'tempsFoulee',  'Stride time',     'SI';
    'distFoulee',   'Stride length',   'SI';
    'DoubleSupport','Double support',  'SI'
};

domains.Variability = {
    'tempsFoulee',  'Stride time',     'CV';
    'distFoulee',   'Stride length',   'CV';
    'DoubleSupport','Double support',  'CV'
};

domainNames = fieldnames(domains);
nDomains    = numel(domainNames);
nGroups     = numel(groups);

% === 2. Compter le nombre total de variables ===
nVarsTotal = 0;
for d = 1:nDomains
    nVarsTotal = nVarsTotal + size(domains.(domainNames{d}), 1);
end

group_means    = zeros(nGroups, nVarsTotal);
group_data     = cell(nGroups, nVarsTotal);
labels         = cell(1, nVarsTotal);
domain_indices = cell(1, nDomains);
stat_type      = cell(1, nVarsTotal);

varIdx = 1;
for d = 1:nDomains
    dName = domainNames{d};
    dVars = domains.(dName);
    nVars = size(dVars, 1);

    domain_indices{d} = varIdx:(varIdx + nVars - 1);

    for v = 1:nVars
        varTech    = dVars{v,1};
        varDisplay = dVars{v,2};
        varType    = dVars{v,3};   % 'Mean', 'CV', 'SI'

        % reconstruire le vrai nom (ex : 'Mean_MoS AP HS (%L0)')
        if isKey(nameMap, varTech)
            baseName = nameMap(varTech);
        else
            baseName = varTech;
        end
       
        varFullName = [varType '_' baseName];
        
% === DEBUG ===
if contains(varFullName, 'MoS AP HS')
    fprintf('🔍 Recherche de : "%s"\n', varFullName);
    if isfield(SpatioTemporalDATA, groups{1}) && ...
       isfield(SpatioTemporalDATA.(groups{1}), condition)
        T = SpatioTemporalDATA.(groups{1}).(condition);
        fprintf('   Variables disponibles contenant "MoS" :\n');
        mos_vars = T.Properties.VariableNames(contains(T.Properties.VariableNames, 'MoS'));
        for i = 1:length(mos_vars)
            fprintf('     - "%s"\n', mos_vars{i});
        end
    end
end
        if strcmp(varType, 'Mean')
            labels{varIdx} = varDisplay;  % pas de mention du type
        else
            labels{varIdx} = sprintf('%s (%s)', varDisplay, varType);
        end

        % collecter pour chaque groupe
        for g = 1:nGroups
            gName = groups{g};

            if isfield(SpatioTemporalDATA, gName) && ...
               isfield(SpatioTemporalDATA.(gName), condition)

                T = SpatioTemporalDATA.(gName).(condition);

                if any(strcmp(T.Properties.VariableNames, varFullName))
                    values = T.(varFullName);
                    group_data{g, varIdx} = values;

                    if numel(values) >= 3 && exist('swtest','file') && swtest(values,0.05)==1
                        group_means(g, varIdx) = median(values, 'omitnan');
                        stat_type{varIdx}      = 'md';
                    else
                        group_means(g, varIdx) = mean(values, 'omitnan');
                        stat_type{varIdx}      = 'μ';
                    end
                else
                    stat_type{varIdx} = ' ';
                end
            end
        end

        varIdx = varIdx + 1;
    end
end

% === 3. Normalisation par variable ===
minVals = zeros(1, nVarsTotal);
maxVals = zeros(1, nVarsTotal);
for v = 1:nVarsTotal
    col_vals = [];
    for g = 1:nGroups
        if ~isempty(group_data{g,v})
            col_vals = [col_vals; group_data{g,v}(:)];
        end
    end
    if ~isempty(col_vals)
        minVals(v) = min(col_vals);
        maxVals(v) = max(col_vals);
    else
        minVals(v) = 0;
        maxVals(v) = 1;
    end
end
scale = @(val,j) (val - minVals(j)) ./ max(maxVals(j) - minVals(j), eps);

% === 4. Angles ===
theta = linspace(0, 2*pi, nVarsTotal+1);
theta(end) = [];
theta = [theta, theta(1)];

% === 5. Figure ===
fig = figure('Name', ['5 Domains Inter - ' condition]);
set(fig, 'Position', FIG_POS);
hold on; axis equal;
set(gca, 'XColor','none', 'YColor','none');

domainColors = {
    [0.7 0.85 1.0];   % Pace
    [1.0 0.7 0.7];    % Rhythm
    [0.7 1.0 0.7];    % Stability
    [1.0 0.85 0.6];   % Asymmetry
    [0.9 0.7 1.0]     % Variability
};

% fonds colorés
for d = 1:nDomains
    idx = domain_indices{d};
    col = domainColors{d};
    for i = idx
        if i < numel(theta)
            patch_theta = [theta(i), theta(i+1), 0];
            patch_r     = [1 1 0];
            [x,y]       = pol2cart(patch_theta, patch_r);
            fill(x, y, col, 'FaceAlpha', 0.45, 'EdgeColor','none', 'HandleVisibility','off');
        end
    end
end

% axes + labels
for i = 1:nVarsTotal
    plot([0 cos(theta(i))], [0 sin(theta(i))], 'k--', 'HandleVisibility','off');

    ang = rad2deg(theta(i));
    x   = LABEL_RADIUS_FACT * cos(theta(i));
    y   = LABEL_RADIUS_FACT * sin(theta(i));
    align = 'left';
    if ang > 90 && ang < 270
        align = 'right';
    end

    label_with_stat = sprintf('%s (%s)', labels{i}, stat_type{i});
    text(x, y, label_with_stat, 'FontSize', 8, ...
        'HorizontalAlignment', align, 'VerticalAlignment','middle', ...
        'Interpreter','none');
end

% couleurs groupes
groupColors = {
    [0 0.6 1];
    [1 0.5 0];
    [0 0.7 0.3];
    [0.5 0 1]
};
groupMarkers = {'^','o','d','s'};

% tracer individus
for g = 1:nGroups
    color = groupColors{g};
    nInd = 0;
    for v = 1:nVarsTotal
        if ~isempty(group_data{g,v})
            nInd = max(nInd, numel(group_data{g,v}));
        end
    end
    for i = 1:nInd
        r_ind = zeros(1, nVarsTotal);
        for v = 1:nVarsTotal
            if ~isempty(group_data{g,v}) && i <= numel(group_data{g,v})
                r_ind(v) = scale(group_data{g,v}(i), v);
            end
        end
        r_ind = [r_ind, r_ind(1)];
        plot(r_ind .* cos(theta), r_ind .* sin(theta), '-', ...
            'Color', [color 0.12], 'LineWidth', 0.6, 'HandleVisibility','off');
    end
end

% tracer moyennes (valeurs en noir)
for g = 1:nGroups
    color  = groupColors{g};
    marker = groupMarkers{g};

    r_norm = arrayfun(@(j) scale(group_means(g,j), j), 1:nVarsTotal);
    r_norm = [r_norm, r_norm(1)];
    x = r_norm .* cos(theta);
    y = r_norm .* sin(theta);

    plot(x, y, '-', 'Color', color, 'LineWidth', 2, ...
        'Marker', marker, 'MarkerSize', 8, 'MarkerFaceColor', color, ...
        'DisplayName', groups{g});

    if SHOW_NUM_VALUES
        for i = 1:nVarsTotal
            if ~isnan(group_means(g,i)) && group_means(g,i) ~= 0
                text(x(i)*1.05, y(i)*1.05, sprintf('%.2f', group_means(g,i)), ...
                    'FontSize', 6.5, 'Color', [0 0 0], ...
                    'HorizontalAlignment','center', 'FontWeight','bold');
            end
        end
    end
end

legend('Location','southoutside', 'Orientation','horizontal', 'FontSize',10);

% ===== blocs de légende à gauche (rapprochés) =====
annotation('textbox', [0.06 0.58 0.14 0.30], 'String', {
    '\bfDomains :', ...
    '\color[rgb]{0.35 0.6 1.0}■ Pace', ...
    '\color[rgb]{1.0 0.3 0.3}■ Rhythm', ...
    '\color[rgb]{0.3 0.8 0.3}■ Stability', ...
    '\color[rgb]{1.0 0.6 0.3}■ Asymmetry', ...
    '\color[rgb]{0.7 0.5 1.0}■ Variability'}, ...
    'FitBoxToText','on', 'BackgroundColor','white', ...
    'EdgeColor','black', 'LineWidth',0.5, 'FontSize',9);

annotation('textbox', [0.06 0.45 0.14 0.10], 'String', {
    '\bfStatistiques :', ...
    '(μ) : Moyenne utilisée', ...
    '(md) : Médiane utilisée'}, ...
    'FitBoxToText','on', 'BackgroundColor','white', ...
    'EdgeColor','black', 'LineWidth',0.5, 'FontSize',9);

annotation('textbox', [0.06 0.23 0.14 0.18], 'String', {
    '\bfGroupes :', ...
    '\color[rgb]{0 0.6 1}— Jeunes Enfants', ...
    '\color[rgb]{1 0.5 0}— Enfants', ...
    '\color[rgb]{0 0.7 0.3}— Adolescents', ...
    '\color[rgb]{0.5 0 1}— Adultes'}, ...
    'FitBoxToText','on', 'BackgroundColor','white', ...
    'EdgeColor','black', 'LineWidth',0.5, 'FontSize',9);

annotation('textbox', [0 0.935 1 0.05], ...
    'String', ['\bf5 Gait Domains - ' condition ' - Inter-group comparison'], ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', ...
    'EdgeColor','none', 'FontSize',14, 'FontWeight','bold');

axis([-1.25 1.25 -1.20 1.30])

end