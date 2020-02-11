% This script runs a comparison between predicted and captured radiance of
% a fluorescent patch.
%
% Copyright, Henryk Blasinski

close all;
clear all;
clc;

ieInit;
wave = 380:4:1000;
nWaves = length(wave);
targetName = 'Target6';

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
saveDir = [];
saveDir = fullfile('~','Desktop','Test images');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% Figure display parameters
fs = 10;
lw = 2;
ms = 7;
sz2 = [1 1 2 2];
sr = 1:5:length(wave);

leg = {'Measured','Predicted'};


% Load illuminant spd's.
illSubset = [];


fName = fullfile(fiToolboxRootPath,'data','Broadband9800K');
illSubset = [illSubset, ieReadSpectra(fName,wave)]; 
% fName = fullfile(fiToolboxRootPath,'data','Broadband6500K');
% illSubset = [illSubset, ieReadSpectra(fName,wave)];
% fName = fullfile(fiToolboxRootPath,'data','Broadband2000K');
% illSubset = [illSubset, ieReadSpectra(fName,wave)];
fName = fullfile(fiToolboxRootPath,'data','NarrowbandBlue');
illSubset = [illSubset, ieReadSpectra(fName,wave)];


fName = fullfile(fiToolboxRootPath,'data','Target6','LEDCubeRadiance','illuminants.mat');
ledCubeLights = ieReadSpectra(fName,wave);


%{
cvx_begin
    variable illMixtureWeights(size(ledCubeLights,2),size(illSubset,2))
    minimize norm(ledCubeLights*illMixtureWeights - illSubset,2)
    subject to
        illMixtureWeights >= 0
cvx_end

figure; plot(illSubset)
hold on;
plot(ledCubeLights*illMixtureWeights ,'--');
%}
% ids = [2, 4, 8, 10, 11];
% ids = [4, 10];
% illSubset = ledCubeLights(:,ids);


dataDir = fullfile(fiToolboxRootPath,'results','experiments');
nCols = 128;
nRows = 103; 

% Load scene reflectance and fluorescence
reflArray = zeros(nRows,nCols,156);
for cc=1:nCols
    
    fName = fullfile(dataDir,sprintf('multiFl_%s_col_%i.mat',targetName,cc));
    data = load(fName);
    
    reflArray(:,cc,:) = data.reflEst';
    
    for rr=1:nRows
        deltaL = data.wave(2) - data.wave(1);
        dm = data.dMatEst{rr}/deltaL;
        flArray(rr,cc) = fluorophoreCreate('type','fromdonaldsonmatrix',...
                                         'DonaldsonMatrix',dm,...
                                         'wave',data.wave);
    end
    
end    


%% Create a scene with the specified reflectance properties and a particular
% illuminant
sceneReTemplate = sceneReflectanceArray(reflArray,1,wave);
sceneFlTemplate = fluorescentSceneCreate('type','fromfluorophore','fluorophore',flArray);


fName = fullfile(fiToolboxRootPath,'data','macbethChart');
refl = ieReadSpectra(fName,wave);
sceneMacbeth = sceneReflectanceArray(reshape(refl',[4 6 nWaves]),10,wave);

% Camera filter transmissivities
fName = fullfile(fiToolboxRootPath,'data','CanonG7x');
cfa = ieReadSpectra(fName,wave);


%% 
delta = 10;

for i=1:size(illSubset,2)      
        
    sceneRe = sceneAdjustIlluminant(sceneReTemplate,illSubset(:,i),0);
    sceneRe = sceneSet(sceneRe,'name',sprintf('refl - ill%i',i));
    ieAddObject(sceneRe);
    
    % Add fluorescence
    sceneReFl = fiSceneAddFluorescence(sceneRe,sceneFlTemplate);
    sceneReFl = sceneSet(sceneReFl,'name',sprintf('refl+fl - ill%i',i));
    ieAddObject(sceneReFl);
    
    
    sceneWindow;


    %% Green
    gYX = [115, 44];

    
    fName = fullfile(fiToolboxRootPath,'data','Target6','LEDCubeRadiance','GreenRadiance.mat');
    measEm = ieReadSpectra(fName,wave);
    measEm = Energy2Quanta(wave,measEm(:,ids(i)));
    measEm = measEm/10^15;
    measEm = measEm(1:156);
    
    estEm = vcGetROIData(sceneReFl,[gYX, delta, delta],'photons');
    estEm = mean(estEm,1);
    estEm = estEm(:);

    % Remove the illuminant from the plot
    estEm(illSubset(:,i) > 0.01*max(illSubset(:,i))) = 0;
    measEm(illSubset(:,i) > 0.01*max(illSubset(:,i))) = 0;
    
    estEm = estEm/max(estEm);
    measEm = measEm/max(measEm);

    figure; 
    hold on; grid on; box on;
    plot(wave,measEm,'r','lineWidth',lw);
    plot(wave,estEm,'g--','lineWidth',lw);
    plot(wave(sr),estEm(sr),'go','markerSize',ms,'lineWidth',lw*0.5);
    plot(wave,illSubset(:,i)/max(illSubset(:,i)),'b');
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);
    
    xlabel('Wavelength, nm','fontsize',fs);
    set(gca,'yTick',0:0.2:1);
    set(gca,'fontsize',fs-2);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPosition',sz2);
    
    if ~isempty(saveDir)
        fName = fullfile(saveDir,sprintf('radiance_green.eps'));
        print('-depsc',fName);
    end


    %% Red
    rYX = [41, 70];

    fName = fullfile(fiToolboxRootPath,'data','Target6','LEDCubeRadiance','RedRadiance.mat');
    measEm = ieReadSpectra(fName,wave);
    measEm = measEm/10^15;
    measEm = measEm(1:156);


    estEm = vcGetROIData(sceneReFl,[rYX delta delta],'photons');
    estEm = mean(estEm,1);
    estEm = estEm(:);

    % Remove the illuminant from the plot
    estEm(illSubset(:,i) > 0.01*max(illSubset(:,i))) = 0;
    measEm(illSubset(:,i) > 0.01*max(illSubset(:,i))) = 0;
    
    estEm = estEm/max(estEm);
    measEm = measEm/max(measEm);

    figure; 
    hold on; grid on; box on;
    plot(wave,measEm,'r','lineWidth',lw);
    plot(wave,estEm,'g--','lineWidth',lw);
    plot(wave(sr),estEm(sr),'go','markerSize',ms,'lineWidth',lw*0.5);
        plot(wave,illSubset(:,i)/max(illSubset(:,i)),'b');

    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);
    
    xlabel('Wavelength, nm','fontsize',fs);
    set(gca,'yTick',0:0.2:1);
    set(gca,'fontsize',fs-2);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPosition',sz2);
    
    if ~isempty(saveDir)
        fName = fullfile(saveDir,sprintf('radiance_red.eps'));
        print('-depsc',fName);
    end


    %% Mix
    mYX = [96, 89];
    
    fName = fullfile(fiToolboxRootPath,'data','Target6','LEDCubeRadiance','MixRadiance.mat');
    measEm = ieReadSpectra(fName,wave);
    measEm = measEm/10^15;
    measEm = measEm(1:156);


    estEm = vcGetROIData(sceneReFl,[mYX delta delta],'photons');
    estEm = mean(estEm,1);
    estEm = estEm(:);

    % Remove the illuminant from the plot
    estEm(illSubset(:,i) > 0.01*max(illSubset(:,i))) = 0;
    measEm(illSubset(:,i) > 0.01*max(illSubset(:,i))) = 0;
    
    estEm = estEm/max(estEm);
    measEm = measEm/max(measEm);
    
    figure; 
    hold on; grid on; box on;
    plot(wave,measEm,'r','lineWidth',lw);
    plot(wave,estEm,'g--','lineWidth',lw);
    plot(wave(sr),estEm(sr),'go','markerSize',ms,'lineWidth',lw*0.5);
    plot(wave,illSubset(:,i)/max(illSubset(:,i)),'b');

    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);
    
    xlabel('Wavelength, nm','fontsize',fs);
    set(gca,'yTick',0:0.2:1);
    set(gca,'fontsize',fs-2);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPosition',sz2);
    
    if ~isempty(saveDir)
        fName = fullfile(saveDir,sprintf('radiance_mix.eps'));
        print('-depsc',fName);
    end

end