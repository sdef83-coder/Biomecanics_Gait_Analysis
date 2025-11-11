%% Voir si problème chez un participant
clc, clear, close all;

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\gaitAnalysisGUI\result\matfiles\ALL')

Participant = {'CTL_24'}; % Exemple 1 seul participant
Condition = {'Plat','Medium','High'};

for iP = 1:length(Participant)
        for iC = 1:length(Condition)
            participant_folder = fullfile('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\gaitAnalysisGUI\result\Fig');
            file = [Participant{iP} '_' Condition{iC} '.mat']; % Charger le fichier pour l'essai actuel
            x = load(file)
            DATA.(Participant{iP}).(Condition{iC}) = x;
        end
        % SECTION GENOU
figure
subplot(3,2,1)
plot(DATA.(Participant{iP}).High.c.data.Left.angle.LKneeAngles  (:,1), 'Color', 'r')
hold on
plot(DATA.(Participant{iP}).Medium.c.data.Left.angle.LKneeAngles(:,1), 'Color', [1, 0.647, 0])
hold on
plot(DATA.(Participant{iP}).Plat.c.data.Left.angle.LKneeAngles(:,1), 'Color', 'g')
grid on
title(sprintf('Left Knee angle %s', Participant{iP}))

subplot(3,2,2)
plot(DATA.(Participant{iP}).High.c.data.Right.angle.RKneeAngles(:,1), 'Color', 'r')
hold on
plot(DATA.(Participant{iP}).Medium.c.data.Right.angle.RKneeAngles(:,1), 'Color', [1, 0.647, 0])
hold on
plot(DATA.(Participant{iP}).Plat.c.data.Right.angle.RKneeAngles(:,1), 'Color', 'g')
grid on
title(sprintf('Right Knee angle %s', Participant{iP}))

% SECTION CHEVILLE
subplot(3,2,3)
plot(DATA.(Participant{iP}).High.c.data.Left.angle.LAnkleAngles  (:,1), 'Color', 'r')
hold on
plot(DATA.(Participant{iP}).Medium.c.data.Left.angle.LAnkleAngles(:,1), 'Color', [1, 0.647, 0])
hold on
plot(DATA.(Participant{iP}).Plat.c.data.Left.angle.LAnkleAngles(:,1), 'Color', 'g')
grid on
title(sprintf('Left Ankle angle %s', Participant{iP}))

subplot(3,2,4)
plot(DATA.(Participant{iP}).High.c.data.Right.angle.RAnkleAngles(:,1), 'Color', 'r')
hold on
plot(DATA.(Participant{iP}).Medium.c.data.Right.angle.RAnkleAngles(:,1), 'Color', [1, 0.647, 0])
hold on
plot(DATA.(Participant{iP}).Plat.c.data.Right.angle.RAnkleAngles(:,1), 'Color', 'g')
grid on
title(sprintf('Right Ankle angle %s', Participant{iP}))


% SECTION HANCHE
subplot(3,2,5)
plot(DATA.(Participant{iP}).High.c.data.Left.angle.LHipAngles(:,1), 'Color', 'r')
hold on
plot(DATA.(Participant{iP}).Medium.c.data.Left.angle.LHipAngles(:,1), 'Color', [1, 0.647, 0])
hold on
plot(DATA.(Participant{iP}).Plat.c.data.Left.angle.LHipAngles(:,1), 'Color', 'g')
grid on
title(sprintf('Left Hip angle %s', Participant{iP}))


subplot(3,2,6)
plot(DATA.(Participant{iP}).High.c.data.Right.angle.RHipAngles(:,1), 'r')
hold on
plot(DATA.(Participant{iP}).Medium.c.data.Right.angle.RHipAngles(:,1), 'Color', [1, 0.647, 0])
hold on
plot(DATA.(Participant{iP}).Plat.c.data.Right.angle.RHipAngles(:,1), 'g')
grid on
title(sprintf('Right Hip angle %s', Participant{iP}))

save_path = fullfile(participant_folder, sprintf('%s_cinematique.png', Participant{iP}));
saveas(gcf, save_path);  % Enregistrer la figure

end