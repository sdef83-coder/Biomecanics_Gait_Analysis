clc, clear, close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\Data\enfants\')

addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\Scripts\btk'));    
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\gaitAnalysisGUI\functions'));

file = ['CTL_09' ...
    '_High_02.c3d'];
acq = btkReadAcquisition(file);

c3d = struct();
c3d.btk = acq;

out = c3dMarkers(c3d);   
%%
data = btkReadAcquisition('CTL_01_High_01.c3d');
markers = btkGetMarkers(data);
plot(markers.RHEE(:,1)); hold on; plot(markers.RHEE(:,2));
legend('X','Y'); title('Déplacement du talon droit');
