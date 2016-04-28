close all;
clear all;
clc;

inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

alpha = 0.01;
beta = 0.01;

maxIter = 1000;

% Create basis function sets
nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = createBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = createBasisSet('excitation','wave',wave','n',nExBasis);
[emBasis, emScore] = createBasisSet('emission','wave',wave','n',nEmBasis);


% Load the light spectra (in photons)
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
illuminant = Energy2Quanta(wave,illuminant);
nChannels = size(illuminant,2);

% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);
       
%% Load simulation data 


nSamples = size(measVals,3);
cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

[ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
fiRecReflAndFl( measVals, camera, cameraGain*deltaL, cameraOffset, illuminant, reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter',maxIter,...
'reflRef',reflRef,'exRef',exRef,'emRef',emRef,'eps',0);

fName = fullfile(fiToolboxRootPath,'results','evaluation','conv_Fl.mat');
save(fName,'reflEst','emEst','exEst','reflValsEst','flValsEst','hist',...
           'maxIter','alpha','beta','exRef','emRef','reflRef','inFName');

       
       
%% Plot the results

% Predicted vs. simulated pixel intensities

figure;
hold all; grid on; box on;
plot(measValsEst(:),measVals(:),'.');
xlabel('Model predicted pixel value');
ylabel('ISET pixel value');

% Convergence
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot([hist{sampleID}.objValsReEm, hist{sampleID}.objValsReEx],'LineWidth',2);
    

end
end

% Estimated vs. ground truth reflectance
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot(wave,reflEst(:,sampleID),'g','LineWidth',2);
    plot(wave,reflRef(:,sampleID),'b--','LineWidth',2);
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);

    rmse = sqrt(mean((reflEst(:,sampleID) - reflRef(:,sampleID)).^2));
    title(sprintf('RMSE %.2f',rmse));

end
end

% Estimated vs. ground truth excitation
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot(wave,exEst(:,sampleID),'g','LineWidth',2);
    plot(wave,exRef(:,sampleID),'b--','LineWidth',2);
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);

    rmse = sqrt(mean((exEst(:,sampleID) - exRef(:,sampleID)).^2));
    title(sprintf('RMSE %.2f',rmse));

end
end

% Estimated vs. ground truth emission
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot(wave,emEst(:,sampleID),'m','LineWidth',2);
    plot(wave,emRef(:,sampleID),'b--','LineWidth',2);
    xlim([min(wave) max(wave)]);

    rmse = sqrt(mean((emEst(:,sampleID) - emRef(:,sampleID)).^2));
    title(sprintf('RMSE %.2f',rmse));

end
end




