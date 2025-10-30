function fig = radarGaitPlot_5Domains_Intra(SpatioTemporalDATA, group, conditions, condColors)

LABEL_RADIUS_FACT = 1.05;
SHOW_NUM_VALUES   = true;
FIG_POS           = [120 100 1450 900];

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
% MoS normalisées
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
    'tempsFoulee',  'Stride time',   'SI';
    'distFoulee',   'Stride length', 'SI';
    'DoubleSupport',  'Double support (%)',  'SI'
};

domains.Variability = {
    'tempsFoulee',  'Stride time',     'CV';
    'distFoulee',   'Stride length',   'CV';
    'DoubleSupport','Double support',  'CV'
};

domainNames = fieldnames(domains);
nDomains = numel(domainNames);
nConds   = numel(conditions);

% compter vars
nVarsTotal = 0;
for d = 1:nDomains
    nVarsTotal = nVarsTotal + size(domains.(domainNames{d}),1);
end

cond_means     = zeros(nConds, nVarsTotal);
cond_data      = cell(nConds, nVarsTotal);
labels         = cell(1, nVarsTotal);
domain_indices = cell(1, nDomains);
stat_type      = cell(1, nVarsTotal);

varIdx = 1;
for d = 1:nDomains
    dVars = domains.(domainNames{d});
    nVars = size(dVars,1);

    domain_indices{d} = varIdx:(varIdx+nVars-1);

    for v = 1:nVars
        varTech    = dVars{v,1};
        varDisplay = dVars{v,2};
        varType    = dVars{v,3};

        if isKey(nameMap, varTech)
            baseName = nameMap(varTech);
        else
            baseName = varTech;
        end

        varFullName = [varType '_' baseName];
        if strcmp(varType, 'Mean')
            labels{varIdx} = varDisplay;  % pas de mention du type
        else
            labels{varIdx} = sprintf('%s (%s)', varDisplay, varType);
        end

        for c = 1:nConds
            cond = conditions{c};
            if isfield(SpatioTemporalDATA, group) && ...
               isfield(SpatioTemporalDATA.(group), cond)

                T = SpatioTemporalDATA.(group).(cond);

                if any(strcmp(T.Properties.VariableNames, varFullName))
                    values = T.(varFullName);
                    cond_data{c,varIdx} = values;

                    if numel(values) >= 3 && exist('swtest','file') && swtest(values,0.05)==1
                        cond_means(c,varIdx) = median(values,'omitnan');
                        stat_type{varIdx}    = 'md';
                    else
                        cond_means(c,varIdx) = mean(values,'omitnan');
                        stat_type{varIdx}    = 'μ';
                    end
                else
                    stat_type{varIdx} = ' ';
                end
            end
        end

        varIdx = varIdx + 1;
    end
end

% normalisation
minVals = zeros(1, nVarsTotal);
maxVals = zeros(1, nVarsTotal);
for v = 1:nVarsTotal
    col_vals = [];
    for c = 1:nConds
        if ~isempty(cond_data{c,v})
            col_vals = [col_vals; cond_data{c,v}(:)];
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
scale = @(val,j) (val - minVals(j)) ./ max(maxVals(j)-minVals(j), eps);

% angles
theta = linspace(0, 2*pi, nVarsTotal+1);
theta(end) = [];
theta = [theta, theta(1)];

fig = figure('Name', ['5 Domains Intra - ' group]);
set(fig, 'Position', FIG_POS);
hold on; axis equal;
set(gca,'XColor','none','YColor','none');

domainColors = {
    [0.7 0.85 1.0];
    [1.0 0.7 0.7];
    [0.7 1.0 0.7];
    [1.0 0.85 0.6];
    [0.9 0.7 1.0]
};

% fonds
for d = 1:nDomains
    idx = domain_indices{d};
    col = domainColors{d};
    for i = idx
        if i < numel(theta)
            patch_theta = [theta(i), theta(i+1), 0];
            patch_r     = [1 1 0];
            [x,y]       = pol2cart(patch_theta, patch_r);
            fill(x,y, col, 'FaceAlpha',0.45, 'EdgeColor','none', 'HandleVisibility','off');
        end
    end
end

% axes + labels
for i = 1:nVarsTotal
    plot([0 cos(theta(i))], [0 sin(theta(i))], 'k--', 'HandleVisibility','off');

    ang = rad2deg(theta(i));
    x   = LABEL_RADIUS_FACT*cos(theta(i));
    y   = LABEL_RADIUS_FACT*sin(theta(i));
    align = 'left';
    if ang > 90 && ang < 270, align = 'right'; end

    label_with_stat = sprintf('%s (%s)', labels{i}, stat_type{i});
    text(x,y, label_with_stat, 'FontSize',8, ...
        'HorizontalAlignment',align, 'VerticalAlignment','middle', ...
        'Interpreter','none');
end

% tracer données individuelles
for c = 1:nConds
    color = condColors{c};
    nInd = 0;
    for v = 1:nVarsTotal
        if ~isempty(cond_data{c,v})
            nInd = max(nInd, numel(cond_data{c,v}));
        end
    end
    for i = 1:nInd
        r_ind = zeros(1, nVarsTotal);
        for v = 1:nVarsTotal
            if ~isempty(cond_data{c,v}) && i <= numel(cond_data{c,v})
                r_ind(v) = scale(cond_data{c,v}(i), v);
            end
        end
        r_ind = [r_ind, r_ind(1)];
        plot(r_ind.*cos(theta), r_ind.*sin(theta), '-', ...
            'Color', [color 0.15], 'LineWidth',0.6, 'HandleVisibility','off');
    end
end

% tracer moyennes
for c = 1:nConds
    color = condColors{c};
    r_norm = arrayfun(@(j) scale(cond_means(c,j), j), 1:nVarsTotal);
    r_norm = [r_norm, r_norm(1)];
    x = r_norm .* cos(theta);
    y = r_norm .* sin(theta);

    plot(x,y, '-', 'Color',color, 'LineWidth',2, ...
        'Marker','o', 'MarkerSize',8, 'MarkerFaceColor',color, ...
        'DisplayName', conditions{c});

    if SHOW_NUM_VALUES
        for i = 1:nVarsTotal
            if ~isnan(cond_means(c,i)) && cond_means(c,i)~=0
                text(x(i)*1.05, y(i)*1.05, sprintf('%.2f', cond_means(c,i)), ...
                    'FontSize',6.5, 'Color', [0 0 0], ...
                    'HorizontalAlignment','center', ...
                    'FontWeight','bold');
            end
        end
    end
end

legend('Location','southoutside', 'Orientation','horizontal', 'FontSize',10);

% ===== blocs de légende à gauche =====
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
    '\bfConditions :', ...
    '\color[rgb]{0.13 0.47 0.85}● Plat', ...
    '\color[rgb]{0 0.6 0}● Medium', ...
    '\color[rgb]{1 0 0}● High'}, ...
    'FitBoxToText','on', 'BackgroundColor','white', ...
    'EdgeColor','black', 'LineWidth',0.5, 'FontSize',9);

annotation('textbox', [0 0.935 1 0.05], ...
    'String', ['\bf5 Gait Domains - ' group ' - Intra-group comparison'], ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', ...
    'EdgeColor','none', 'FontSize',14, 'FontWeight','bold');

axis([-1.25 1.25 -1.20 1.30])

end