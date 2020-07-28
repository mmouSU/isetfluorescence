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

% NADH
ieNewGraphWin;
NADHEmission = ieReadSpectra('NADPH_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,NADHEmission,'b','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
NADH = fiReadFluorophore('NADH.mat','wave',theseWaves); 
plot(theseWaves,NADH.emission,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
NADH = fiReadFluorophore('NADH.mat','wave',theseWaves); 
plot(theseWaves,NADH.emission,'r','LineWidth',3); hold on;
legend('Palero Et al 2007', 'Monici 2005', 'DaCosta et al 2003');
title('NADH Emission');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');

% collagen
ieNewGraphWin;
collagenEmission = ieReadSpectra('collagen_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,collagenEmission,'b','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
collagen = fiReadFluorophore('Collagen.mat','wave',theseWaves); 
plot(theseWaves,collagen.emission,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
collagen = fiReadFluorophore('collagen1.mat','wave',theseWaves); 
plot(theseWaves,collagen.emission,'r','LineWidth',3); hold on;
legend('Palero Et al 2007', 'Monici 2005', 'DaCosta et al 2003');
title('Collagen Emission');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');

% elastin 
ieNewGraphWin;
elastinEmission = ieReadSpectra('elastin_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,elastinEmission,'b','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
elastin = fiReadFluorophore('Elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.emission,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
elastin = fiReadFluorophore('elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.emission,'r','LineWidth',3); hold on;
legend('Palero Et al 2007', 'Monici 2005', 'DaCosta et al 2003');
title('Elastin Emission');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');

% FAD
ieNewGraphWin;
FADEmission = ieReadSpectra('FAD_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,FADEmission,'b','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.emission,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.emission,'r','LineWidth',3); hold on;
legend('Palero Et al 2007', 'Monici 2005', 'DaCosta et al 2003');
title('FAD Emission');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');

% Porphyrins
ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
Porphyrins = fiReadFluorophore('Porphyrins.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.emission,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
Porphyrins= fiReadFluorophore('protoporphyrin.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.emission,'r','LineWidth',3); hold on;
legend('Monici 2005', 'DaCosta et al 2003');
title('Porphyrin Emission');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');


%% Plot Excitation

ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
NADH = fiReadFluorophore('NADH.mat','wave',theseWaves); 
plot(theseWaves,NADH.excitation,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
NADH = fiReadFluorophore('NADH.mat','wave',theseWaves); 
plot(theseWaves,NADH.excitation,'r','LineWidth',3); hold on;
legend('Monici 2005', 'DaCosta et al 2003');
title('NADH Excitation');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');
x = [385,385]; y = [0,1]; plot(x,y,'k--','LineWidth',3);

ieNewGraphWin;
% Deal et al 2018 does not have a "real" maximum
% chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
% collagenExcitation = ieReadSpectra('collagen_excitation_DealEtAl2018.mat',theseWaves); 
% plot(theseWaves,collagenExcitation,'b'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
collagen = fiReadFluorophore('Collagen.mat','wave',theseWaves); 
plot(theseWaves,collagen.excitation,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
collagen = fiReadFluorophore('collagen1.mat','wave',theseWaves); 
plot(theseWaves,collagen.excitation,'r','LineWidth',3); hold on;
legend('Monici 2005', 'DaCosta et al 2003');
title('Collagen Excitation');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');

ieNewGraphWin;
% Deal et al 2018 does not have a "real" maximum
% chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
% elastinExcitation = ieReadSpectra('elastin_excitation_DealEtAl2018.mat',theseWaves); 
% plot(theseWaves,collagenExcitation,'b'); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
elastin = fiReadFluorophore('Elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.excitation,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
elastin = fiReadFluorophore('elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.excitation,'r','LineWidth',3); hold on;
legend('Monici 2005', 'DaCosta et al 2003');
title('Elastin Excitation');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');

ieNewGraphWin;
FADExcitation = ieReadSpectra('FAD_excitation_DealEtAl2018.mat',theseWaves); 
plot(theseWaves,FADExcitation,'b','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.excitation,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.excitation,'r','LineWidth',3); hold on;
legend('Deal et al 2019', 'Monici 2005', 'DaCosta et al 2003');
title('FAD Excitation');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');

ieNewGraphWin;
PorphyrinsExcitation = ieReadSpectra('porphyrins_excitation_DealEtAl2018.mat',theseWaves); 
plot(theseWaves,PorphyrinsExcitation,'b','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
Porphyrins = fiReadFluorophore('Porphyrins.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.excitation,'g','LineWidth',3); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
Porphyrins = fiReadFluorophore('protoporphyrin.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.excitation,'r','LineWidth',3); hold on;
legend('Deal et al 2019', 'Monici 2005', 'DaCosta et al 2003');
title('Porphyrins Excitation');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');


%%
fullFileName = which('oeCameraQE.mat');
[data,filterNames] = ieReadColorFilter(wave,fullFileName); 
ieNewGraphWin;
plot(wave,data);