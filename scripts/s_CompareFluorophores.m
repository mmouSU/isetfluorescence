% s_CompareFluorophores.m

% Read in the excitation and emission spectra for different fluorophores from different
% sources and compare
% The best data will have the entire spectra mapped out ... otherwise the
% max which is set to 1 will not be the actual maximum excitation or
% emission
% For this reason, we discard any data for which the entire spectra is not
% measured (e.g. Deal et al 2018)

% They are close enough
theseWaves = 300:5:700;

% Plot Emission Spectra from 3 different sources
% collagen
ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
collagenEmission = ieReadSpectra('collagen_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,collagenEmission,'b'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
collagen = fiReadFluorophore('Collagen.mat','wave',theseWaves); 
plot(theseWaves,collagen.emission,'g'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
collagen = fiReadFluorophore('collagen1.mat','wave',theseWaves); 
plot(theseWaves,collagen.emission,'r','LineWidth',3); hold on;
legend('Palero Et al 2007', 'Monici', 'DaCosta et al 2003');
title('Collagen Emission');
ax = gca;
ax.FontSize=16;

% elastin 
ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
elastinEmission = ieReadSpectra('elastin_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,elastinEmission,'b'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
elastin = fiReadFluorophore('Elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.emission,'g'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
elastin = fiReadFluorophore('elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.emission,'r','LineWidth',3); hold on;
legend('Palero Et al 2007', 'Monici', 'DaCosta et al 2003');
title('Elastin Emission');
ax = gca;
ax.FontSize=16;

% FAD
ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
FADEmission = ieReadSpectra('FAD_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,FADEmission,'b'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.emission,'g'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.emission,'r','LineWidth',3); hold on;
legend('Palero Et al 2007', 'Monici', 'DaCosta et al 2003');
title('FAD Emission');
ax = gca;
ax.FontSize=16;

% Porphyrins
ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
Porphyrins = fiReadFluorophore('Porphyrins.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.emission,'g'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
Porphyrins= fiReadFluorophore('protoporphyrin.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.emission,'r','LineWidth',3); hold on;
legend('Monici', 'DaCosta et al 2003');
title('Porphyrin Emission');
ax = gca;
ax.FontSize=16;


%% Plot Excitation

ieNewGraphWin;
% Deal et al 2018 does not have a "real" maximum
% chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
% collagenExcitation = ieReadSpectra('collagen_excitation_DealEtAl2018.mat',theseWaves); 
% plot(theseWaves,collagenExcitation,'b'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
collagen = fiReadFluorophore('Collagen.mat','wave',theseWaves); 
plot(theseWaves,collagen.excitation,'g'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
collagen = fiReadFluorophore('collagen1.mat','wave',theseWaves); 
plot(theseWaves,collagen.excitation,'r','LineWidth',3); hold on;
legend('Monici', 'DaCosta et al 2003');
title('Collagen Excitation');
ax = gca;
ax.FontSize=16;

ieNewGraphWin;
% Deal et al 2018 does not have a "real" maximum
% chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
% elastinExcitation = ieReadSpectra('elastin_excitation_DealEtAl2018.mat',theseWaves); 
% plot(theseWaves,collagenExcitation,'b'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
elastin = fiReadFluorophore('Elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.excitation,'g'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
elastin = fiReadFluorophore('elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.excitation,'r','LineWidth',3); hold on;
legend('Monici', 'DaCosta et al 2003');
title('Elastin Excitation');
ax = gca;
ax.FontSize=16;

ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
FADExcitation = ieReadSpectra('FAD_excitation_DealEtAl2018.mat',theseWaves); 
plot(theseWaves,FADExcitation,'b'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.excitation,'g'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.excitation,'r','LineWidth',3); hold on;
legend('Deal et al 2019', 'Monici', 'DaCosta et al 2003');
title('FAD Excitation');
ax = gca;
ax.FontSize=16;

ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
PorphyrinsExcitation = ieReadSpectra('porphyrins_excitation_DealEtAl2018.mat',theseWaves); 
plot(theseWaves,PorphyrinsExcitation,'b'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
Porphyrins = fiReadFluorophore('Porphyrins.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.excitation,'g'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
Porphyrins = fiReadFluorophore('protoporphyrin.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.excitation,'r','LineWidth',3); hold on;
legend('Deal et al 2019', 'Monici', 'DaCosta et al 2003');
title('Porphyrins Excitation');
ax = gca;
ax.FontSize=16;

