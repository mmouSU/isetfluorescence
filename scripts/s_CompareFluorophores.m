% s_CompareFluorophores.m

% Read in the excitation and emission spectra for different fluorophores from different
% sources and compare

% They are close enough
theseWaves = 300:5:700;

% Plot Emission Spectra from 3 different sources
% collagen
ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
collagenEmission = ieReadSpectra('collagen_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,collagenEmission); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
collagen = fiReadFluorophore('Collagen.mat','wave',theseWaves); 
plot(theseWaves,collagen.emission); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
collagen = fiReadFluorophore('collagen1.mat','wave',theseWaves); 
plot(theseWaves,collagen.emission); hold on;

% elastin 
ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
elastinEmission = ieReadSpectra('elastin_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,elastinEmission); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
elastin = fiReadFluorophore('Elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.emission); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
elastin = fiReadFluorophore('elastin.mat','wave',theseWaves); 
plot(theseWaves,elastin.emission); hold on;

% FAD
ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','OtherSources'));
FADEmission = ieReadSpectra('FAD_emission_PaleroEtAl2007.mat',theseWaves); 
plot(theseWaves,FADEmission); hold on;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.emission); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
FAD = fiReadFluorophore('FAD.mat','wave',theseWaves); 
plot(theseWaves,FAD.emission); hold on;

% Porphyrins
ieNewGraphWin;
chdir(fullfile(fiToolboxRootPath,'data','Monici'));
Porphyrins = fiReadFluorophore('Porphyrins.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.emission); hold on;
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
Porphyrins= fiReadFluorophore('protoporphyrin.mat','wave',theseWaves); 
plot(theseWaves,Porphyrins.emission); hold on;


%% Plot Excitation

% elastinEmission = ieReadSpectra('elastin_emission_PaleroEtAl2007.mat',theseWaves); 
% plot(theseWaves,elastinEmission); 
% FADEmission = ieReadSpectra('FAD_emission_PaleroEtAl2007.mat',theseWaves); 
% plot(theseWaves,FADEmission); 
% NADPHEmission = ieReadSpectra('NADPH_emission_PaleroEtAl2007.mat',theseWaves); 
% plot(theseWaves,NADPHEmission);
% title('Emission')
% legend('Collagen Palero Et Al 2007', 'Elastin Palero Et Al 2007','FAD Palero Et Al 2007','NADPH Palero Et Al 2007');
% 
% 
% 
% ieNewGraphWin;
% collagenExcitation = ieReadSpectra('collagen_excitation_DealEtAl2018.mat',theseWaves); 
% plot(theseWaves,collagenExcitation * .35,'m'); hold on;
% elastinExcitation  = ieReadSpectra('elastin_excitation_DealEtAl2018.mat',theseWaves); 
% plot(theseWaves,elastinExcitation,'r'); 
% FADExcitation  = ieReadSpectra('FAD_excitation_DealEtAl2018.mat',theseWaves); 
% plot(theseWaves,FADExcitation); 
% NADHExcitation  = ieReadSpectra('NADH_excitation_DealEtAl2018.mat',theseWaves); 
% plot(theseWaves,NADHExcitation * .45,'b');
% porphyrinsExcitation = ieReadSpectra('porphyrins_excitation_DealEtAl2018.mat',theseWaves); 
% plot(theseWaves,porphyrinsExcitation);
% title('Excitation')
% legend('Collagen Deal Et al 2018', 'Elastin Deal Et al 2018','FAD Deal Et al 2018','NADH Deal Et al 2018','Porphyrins Deal Et al 2018');

 
 


 
 
 
