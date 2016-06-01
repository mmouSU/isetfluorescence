% Use the CIM algorithm to estimate the reflectance and fluorescence
% properties from simulated data.
%
% Copyright, Henryk Blasinski 2016

close all;
clear variables;
clc;

% Load data
inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

% Wavelength sampling is defined in the loaded data file
% wave = 400:4:1000;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

alpha = 0.1;
beta = 0.1;

% Create basis function sets
nReflBasis = 5;
nEmBasis = 12;

[reflBasis, reflScore] = fiCreateBasisSet('reflectance','wave',wave','n',nReflBasis);
[emBasis, emScore] = fiCreateBasisSet('emission','wave',wave','n',nEmBasis);

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
      
%% Perform estimation

nSamples = size(measVals,3);
cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

[ reflEst, rfCoeffs, emEst, emCoeffs, emWghts, reflValsEst, flValsEst, hist ] = ...
fiRecReflAndEm( measVals, camera, cameraGain*deltaL, cameraOffset, illuminant, reflBasis, emBasis, alpha, beta, 'maxIter',10);

measValsEst = flValsEst + reflValsEst;


%% Plot the results

[err, std] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals - cameraOffset,[nChannels*nFilters,nSamples]), 'absolute');
fprintf('Total pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), 'absolute');
fprintf('Reflected pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), 'absolute');
fprintf('Fluoresced pixel error %.3f, std %.3f\n',err,std);


[err, std] = fiComputeError(reflEst, reflRef, 'absolute');
fprintf('Reflectance error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(emEst, emRef, 'normalized');
fprintf('Emission error (normalized) %.3f, std %.3f\n',err,std);


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
    plot([hist{sampleID}.objValsReEm, hist{sampleID}.objValsReWe],'LineWidth',2);
    

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


% Estimated vs. ground truth normalized emission
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot(wave,emEst(:,sampleID)/max(emEst(:,sampleID)),'m','LineWidth',2);
    plot(wave,emRef(:,sampleID)/max(emRef(:,sampleID)),'b--','LineWidth',2);
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);

    rmse = sqrt(mean((emEst(:,sampleID)/max(emEst(:,sampleID)) - emRef(:,sampleID)/max(emRef(:,sampleID))).^2));
    title(sprintf('RMSE %.2f',rmse));

end
end



