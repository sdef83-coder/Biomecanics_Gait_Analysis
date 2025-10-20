function [stats] = Spatiotempocalc(Left, Right)
% CycleVariability - Calcule moyennes, médianes, SD, CV et SI pour plusieurs variables
%
% Entrée :
%   Left, Right : structures contenant les cycles pour chaque jambe
%
% Sortie :
%   stats : structure contenant les statistiques par jambe + moyennes + médianes + indices de symétrie

    variableNames = {
        'tempsFoulee', ...
        'vitCadencePasParMinute', ...
        'vitFoulee', ...
        'pctToeOffOppose', ...
        'pctSimpleAppuie', ...
        'distFoulee', ...
        'comRangeHeight', ...
        'comRangeML', ...
        'NormStepLength', ...
        'NormCadence', ...
        'NormWalkRatio', ...
        'LargeurPas', ...
        'MoS_AP_Mean', ...
        'MoS_ML_Mean', ...
        'MoS_HS_AP', ...
        'MoS_HS_ML'
    };

    for v = 1:length(variableNames)
        varName = variableNames{v};

        leftVar  = vertcat(Left.(varName));
        rightVar = vertcat(Right.(varName));

        % Conversion mm → cm pour comRangeML
        if strcmp(varName, 'comRangeML')
            leftVar = leftVar / 10;
            rightVar = rightVar / 10;
        end

        if strcmp(varName, 'vitCadencePasParMinute')
            leftVar  = 2 * leftVar;   % conversion foulées/min → pas/min
            rightVar = 2 * rightVar;
        end

        % Moyennes
        stats.([varName '_Mean_Left']) = mean(leftVar, 'omitnan');
        stats.([varName '_Mean_Right']) = mean(rightVar, 'omitnan');
        stats.([varName '_Mean_Mean']) = mean([ ...
            stats.([varName '_Mean_Left']), ...
            stats.([varName '_Mean_Right'])], 'omitnan');

        % Écarts-types
        stats.([varName '_SD_Left']) = std(leftVar, 'omitnan');
        stats.([varName '_SD_Right']) = std(rightVar, 'omitnan');

        % Coefficients de variation (basés sur les moyennes)
        stats.([varName '_CV_Left']) = ...
            (stats.([varName '_SD_Left']) / stats.([varName '_Mean_Left'])) * 100;
        stats.([varName '_CV_Right']) = ...
            (stats.([varName '_SD_Right']) / stats.([varName '_Mean_Right'])) * 100;
        stats.([varName '_CV_Mean']) = mean([ ...
            stats.([varName '_CV_Left']), ...
            stats.([varName '_CV_Right'])], 'omitnan');

        % Symmetry Index (basé sur les moyennes)
        L = stats.([varName '_Mean_Left']);
        R = stats.([varName '_Mean_Right']);
        stats.([varName '_SI']) = (abs(L - R) / (0.5 * (L + R))) * 100;
    end

    % === Calcul du double appui pour chaque cycle ===
    % Gauche
    left_DS = vertcat(Left.pctToeOffOppose) + ...
              (vertcat(Left.pctToeOff) - vertcat(Left.pctContactTalOppose));
    % Droite
    right_DS = vertcat(Right.pctToeOffOppose) + ...
               (vertcat(Right.pctToeOff) - vertcat(Right.pctContactTalOppose));

    % Moyennes
    stats.DoubleSupport_Mean_Left = mean(left_DS, 'omitnan');
    stats.DoubleSupport_Mean_Right = mean(right_DS, 'omitnan');
    stats.DoubleSupport_Mean_Mean = mean([ ...
        stats.DoubleSupport_Mean_Left, ...
        stats.DoubleSupport_Mean_Right], 'omitnan');

    % Écarts-types
    stats.DoubleSupport_SD_Left = std(left_DS, 'omitnan');
    stats.DoubleSupport_SD_Right = std(right_DS, 'omitnan');

    % Coefficients de variation (basés sur les moyennes)
    stats.DoubleSupport_CV_Left = ...
        (stats.DoubleSupport_SD_Left / stats.DoubleSupport_Mean_Left) * 100;
    stats.DoubleSupport_CV_Right = ...
        (stats.DoubleSupport_SD_Right / stats.DoubleSupport_Mean_Right) * 100;
    stats.DoubleSupport_CV_Mean = mean([ ...
        stats.DoubleSupport_CV_Left, stats.DoubleSupport_CV_Right], 'omitnan');

    % Symmetry Index (basé sur les moyennes)
    L = stats.DoubleSupport_Mean_Left;
    R = stats.DoubleSupport_Mean_Right;
    stats.DoubleSupport_SI = (abs(L - R) / (0.5 * (L + R))) * 100;

    % === Calcul du Walk Ratio par jambe (en cm/pas/min) ===
    % Gauche
    leftStrideLength  = vertcat(Left.distPas);
    leftCadence       = 2*vertcat(Left.vitCadencePasParMinute);
    leftWalkRatio     = (leftStrideLength ./ leftCadence) * 100;  % en cm·min⁻¹·pas⁻¹

    % Droite
    rightStrideLength = vertcat(Right.distPas);
    rightCadence      = 2*vertcat(Right.vitCadencePasParMinute);
    rightWalkRatio    = (rightStrideLength ./ rightCadence) * 100;  % en cm·min⁻¹·pas⁻¹

    % Moyennes
    stats.WalkRatio_Mean_Left = mean(leftWalkRatio, 'omitnan');
    stats.WalkRatio_Mean_Right = mean(rightWalkRatio, 'omitnan');
    stats.WalkRatio_Mean_Mean = mean([ ...
        stats.WalkRatio_Mean_Left, ...
        stats.WalkRatio_Mean_Right], 'omitnan');

    % Écarts-types
    stats.WalkRatio_SD_Left = std(leftWalkRatio, 'omitnan');
    stats.WalkRatio_SD_Right = std(rightWalkRatio, 'omitnan');

    % Coefficients de variation (basés sur les moyennes)
    stats.WalkRatio_CV_Left = ...
        (stats.WalkRatio_SD_Left / stats.WalkRatio_Mean_Left) * 100;
    stats.WalkRatio_CV_Right = ...
        (stats.WalkRatio_SD_Right / stats.WalkRatio_Mean_Right) * 100;
    stats.WalkRatio_CV_Mean = mean([ ...
        stats.WalkRatio_CV_Left, stats.WalkRatio_CV_Right], 'omitnan');

    % Symmetry Index (basé sur les moyennes)
    L = stats.WalkRatio_Mean_Left;
    R = stats.WalkRatio_Mean_Right;
    stats.WalkRatio_SI = (abs(L - R) / (0.5 * (L + R))) * 100;
end