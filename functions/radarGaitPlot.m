function fig = radarGaitPlot(SpatioTemporalDATA, condition, param_names)

% === Données
nParams = length(param_names);
N_jeune = length(SpatioTemporalDATA.JeunesEnfants.(condition).(['Mean_' param_names{1}]));
N_enf = length(SpatioTemporalDATA.Enfants.(condition).(['Mean_' param_names{1}]));
N_ado = length(SpatioTemporalDATA.Adolescents.(condition).(['Mean_' param_names{1}]));
N_adu = length(SpatioTemporalDATA.Adultes.(condition).(['Mean_' param_names{1}]));
data_jeunes = zeros(N_jeune, nParams*2);
data_enfants = zeros(N_enf, nParams*2);
data_ado = zeros(N_ado, nParams*2);
data_adultes = zeros(N_adu, nParams*2);
for i = 1:nParams
    base = param_names{i};
    data_jeunes(:,i) = SpatioTemporalDATA.JeunesEnfants.(condition).(['Mean_' base]);
    data_jeunes(:,i+nParams) = SpatioTemporalDATA.JeunesEnfants.(condition).(['CV_' base]);
    data_enfants(:,i) = SpatioTemporalDATA.Enfants.(condition).(['Mean_' base]);
    data_enfants(:,i+nParams) = SpatioTemporalDATA.Enfants.(condition).(['CV_' base]);
    data_ado(:,i) = SpatioTemporalDATA.Adolescents.(condition).(['Mean_' base]);
    data_ado(:,i+nParams) = SpatioTemporalDATA.Adolescents.(condition).(['CV_' base]);
    data_adultes(:,i) = SpatioTemporalDATA.Adultes.(condition).(['Mean_' base]);
    data_adultes(:,i+nParams) = SpatioTemporalDATA.Adultes.(condition).(['CV_' base]);
end

% === Normalisation
all_data = [data_jeunes; data_enfants; data_ado; data_adultes];
minVals = min(all_data,[],1);
maxVals = max(all_data,[],1);
scale_data = @(val, minV, maxV) (val - minV) / (maxV - minV);

% === Moyenne ou médiane selon normalité
nVar = size(all_data,2);
mean_jeune = zeros(1, nVar);
mean_enf = zeros(1, nVar);
mean_ado = zeros(1, nVar);
mean_adu = zeros(1, nVar);
stat_labels = cell(1, nVar); % (μ) ou (md)

for j = 1:nVar
    % Enfants
    if swtest(data_enfants(:,j), 0.05) == 1
        mean_enf(j) = median(data_enfants(:,j), 'omitnan');
        stat_labels{j} = '(md)';
    else
        mean_enf(j) = mean(data_enfants(:,j), 'omitnan');
        stat_labels{j} = '(μ)';
    end

    % Adultes (on ne change pas le label si déjà médiane)
    if swtest(data_adultes(:,j), 0.05) == 1
        mean_adu(j) = median(data_adultes(:,j), 'omitnan');
        if ~strcmp(stat_labels{j}, '(md)')
            stat_labels{j} = '(md)';
        end
    else
        mean_adu(j) = mean(data_adultes(:,j), 'omitnan');
    end

    % Adolescents
 if swtest(data_ado(:,j), 0.05) == 1
        mean_ado(j) = median(data_ado(:,j), 'omitnan');
        if ~strcmp(stat_labels{j}, '(md)')
            stat_labels{j} = '(md)';
        end
    else
        mean_ado(j) = mean(data_ado(:,j), 'omitnan');
    end

   % Jeunes enfants
      if swtest(data_jeunes(:,j), 0.05) == 1
        mean_jeune(j) = median(data_jeunes(:,j), 'omitnan');
        if ~strcmp(stat_labels{j}, '(md)')
            stat_labels{j} = '(md)';
        end
    else
        mean_jeune(j) = mean(data_jeunes(:,j), 'omitnan');
    end
end

% === Création des labels
labels = cell(1, nVar);
for i = 1:nParams
    base = param_names{i};
    nom = base;
    labels{i} = [nom ' ' stat_labels{i}];
    labels{i+nParams} = [extractBefore(nom, ' (') ' (CV) ' stat_labels{i+nParams}];
end

% === Radar params
theta = linspace(0, 2*pi, nVar+1); theta(end) = []; theta = [theta, theta(1)];
color_jeune = [0 0.6 1]; color_enf = [1, 0.5, 0]; color_ado = [0 0.70 0.30]; color_adu = [0.5, 0, 1];

% === Figure
fig = figure;
hold on; axis equal;
set(gca, 'XColor', 'none', 'YColor', 'none');

% === Zones colorées
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
        fill(x, y, col, 'FaceAlpha', 0.25, 'EdgeColor','none');
    end
end

% === Axes + labels
for i = 1:nVar
    x = [0, cos(theta(i))]; y = [0, sin(theta(i))];
    plot(x, y, 'k--');
    xpos = 1.05 * cos(theta(i)); ypos = 1.05 * sin(theta(i));
    align = 'left'; if rad2deg(theta(i)) > 90 && rad2deg(theta(i)) < 270, align = 'right'; end
    text(xpos, ypos, labels{i}, 'HorizontalAlignment', align, ...
        'FontSize', 8, 'Rotation', 0, 'VerticalAlignment', 'middle');
end

% === Courbes
for i = 1:N_jeune
    r = arrayfun(@(j) scale_data(data_jeunes(i,j), minVals(j), maxVals(j)), 1:nVar);
    r = [r, r(1)];
    plot(r .* cos(theta), r .* sin(theta), '-', 'Color', [color_jeune 0.3]);
end
for i = 1:N_enf
    r = arrayfun(@(j) scale_data(data_enfants(i,j), minVals(j), maxVals(j)), 1:nVar);
    r = [r, r(1)];
    plot(r .* cos(theta), r .* sin(theta), '-', 'Color', [color_enf 0.3]);
end
for i = 1:N_ado
    r = arrayfun(@(j) scale_data(data_ado(i,j), minVals(j), maxVals(j)), 1:nVar);
    r = [r, r(1)];
    plot(r .* cos(theta), r .* sin(theta), '-', 'Color', [color_ado 0.3]);
end
for i = 1:N_adu
    r = arrayfun(@(j) scale_data(data_adultes(i,j), minVals(j), maxVals(j)), 1:nVar);
    r = [r, r(1)];
    plot(r .* cos(theta), r .* sin(theta), '-', 'Color', [color_adu 0.3]);
end

% === Moyenne/médiane jeunes enfants
r_jeune = arrayfun(@(j) scale_data(mean_jeune(j), minVals(j), maxVals(j)), 1:nVar);
r_jeune = [r_jeune, r_jeune(1)];
x_jeune = r_jeune .* cos(theta); y_jeune = r_jeune .* sin(theta);
plot(x_jeune, y_jeune, '-', 'Color', color_jeune, 'LineWidth', 2, 'Marker', '^');

% === Moyenne/médiane enfants
r_enf = arrayfun(@(j) scale_data(mean_enf(j), minVals(j), maxVals(j)), 1:nVar);
r_enf = [r_enf, r_enf(1)];
x_enf = r_enf .* cos(theta); y_enf = r_enf .* sin(theta);
plot(x_enf, y_enf, '-', 'Color', color_enf, 'LineWidth', 2, 'Marker', 'o');

% === Moyenne/médiane adolescents
r_ado = arrayfun(@(j) scale_data(mean_ado(j), minVals(j), maxVals(j)), 1:nVar);
r_ado = [r_ado, r_ado(1)];
x_ado = r_ado .* cos(theta); y_ado = r_ado .* sin(theta);
plot(x_ado, y_ado, '-', 'Color', color_ado, 'LineWidth', 2, 'Marker', 'd');

% === Moyenne/médiane adultes
r_adu = arrayfun(@(j) scale_data(mean_adu(j), minVals(j), maxVals(j)), 1:nVar);
r_adu = [r_adu, r_adu(1)];
x_adu = r_adu .* cos(theta); y_adu = r_adu .* sin(theta);
plot(x_adu, y_adu, '-', 'Color', color_adu, 'LineWidth', 2, 'Marker', 's');

% === Valeurs numériques
for i = 1:nVar
    text(x_jeune(i)*1.05, y_jeune(i)*1.05, sprintf('%.2f', mean_jeune(i)), 'Color', 'k', 'FontSize', 8);
    text(x_enf(i)*1.05, y_enf(i)*1.05, sprintf('%.2f', mean_enf(i)), 'Color', 'k', 'FontSize', 8);
    text(x_ado(i)*1.05, y_ado(i)*1.05, sprintf('%.2f', mean_ado(i)), 'Color', 'k', 'FontSize', 8);
    text(x_adu(i)*1.05, y_adu(i)*1.05, sprintf('%.2f', mean_adu(i)), 'Color', 'k', 'FontSize', 8);
end

% === Légendes
% Légende des groupes
annotation('textbox', [0.01 0.26 0.3 0.25], 'String', {
    '\bfGroupes :', ...
    '\color[rgb]{0 0.6 1}⎯⎯⎯⎯⎯ Jeunes Enfants (individuels)', ...
    '\color[rgb]{1 0.5 0}⎯⎯⎯⎯⎯ Enfants (individuels)', ...
    '\color[rgb]{0 0.70 0.30}⎯⎯⎯⎯⎯ Adolescents (individuels)', ...
    '\color[rgb]{0.5 0 1}⎯⎯⎯⎯⎯ Adultes (individuels)', ...
    '\color[rgb]{0 0.6 1}▲ Moyenne/Médiane Jeunes Enfants', ...
    '\color[rgb]{1 0.5 0}● Moyenne/Médiane Enfants', ...
    '\color[rgb]{0 0.70 0.30}◆ Moyenne/Médiane Adolescents', ...
    '\color[rgb]{0.5 0 1}■ Moyenne/Médiane Adultes'}, ...
    'FitBoxToText','on', ...
    'EdgeColor','none', ...
    'FontSize', 9);

% Légende Statistiques
annotation('textbox', [0.01 0.55 0.3 0.05], 'String', {
    '\bfStatistiques :', ...
    '(μ) : Moyenne utilisée', ...
    '(md) : Médiane utilisée'}, ...
    'FitBoxToText','on', ...
    'EdgeColor','none', ...
    'FontSize', 9);

% Légende Couleurs de fond (zones Pace, Rhythm...)
annotation('textbox', [0.01 0.70 0.2 0.15], 'String', {
    '\color[rgb]{0.35 0.6 1.0} Pace', ...
    '\color[rgb]{1.0 0.3 0.3} Rhythm', ...
    '\color[rgb]{0.3 0.8 0.3} Stability', ...
    '\color[rgb]{0.85 0.85 0} Variability'}, ...
    'FitBoxToText','on', 'BackgroundColor','none', ...
    'EdgeColor','none', 'FontSize', 9);

annotation('textbox', [0 0.955 1 0.05], ...
    'String', ['\bfSpatiotemporal variables - ' condition], ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'EdgeColor', 'none', 'FontSize', 14);

end