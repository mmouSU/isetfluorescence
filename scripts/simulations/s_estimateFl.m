% Use the single fluorophore algorithm to estimate the reflectance and fluorescence
% properties from simulated data.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Load simulation data
inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

% Wavelength sampling is defined in the simulation data file.
deltaL = wave(2) - wave(1);
nWaves = length(wave);

alpha = 0.01;
beta = 0.01;

% Create basis function sets
nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = fiCreateBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = fiCreateBasisSet('excitation','wave',wave','n',nExBasis);
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
       
%% Load simulation data 

nSamples = size(measVals,3);
% cameraGain = repmat(cameraGain,[1 1 nSamples]);
% cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

measValsNoise = max(measVals + 0.01*randn(size(measVals)),0);
        
localCameraGain = repmat(cameraGain,[1 1 nSamples]);
localCameraOffset = repmat(cameraOffset,[1 1 nSamples]);
        
nF = max(max(measValsNoise,[],1),[],2);
localCameraGain = localCameraGain./repmat(nF,[nFilters nChannels 1]);
measValsNoise = measValsNoise./repmat(nF,[nFilters nChannels 1]);


for i=1:3
    switch i
        case 1
            alpha = 0;
            beta = 0.1;
        case 2
            alpha = 0.1;
            beta = 0;
        case 3
            alpha = 0.1;
            beta = 0.1;
    end
    
    [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndFl( measValsNoise, camera, localCameraGain*deltaL, localCameraOffset, illuminant, reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter', 20,...
        'reflRef',reflRef,'exRef',exRef,'emRef',emRef,'eps',0);
    
    
    measValsEst = reflValsEst + flValsEst;
    
    fprintf('====== alpha=%.3f beta=%.3f ======\n',alpha, beta);
    
    fprintf('Total time %f, per sample %f, per iteration %f\n',sum(cellfun(@(x) sum(x.computeTime),hist)),mean(cellfun(@(x) sum(x.computeTime),hist)),mean(cellfun(@(x) mean(x.computeTime),hist)));
    
    [err, std] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals - cameraOffset,[nChannels*nFilters,nSamples]), 'absolute');
    fprintf('Total pixel error %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), 'absolute');
    fprintf('Reflected pixel error %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), 'absolute');
    fprintf('Fluoresced pixel error %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(reflEst, reflRef, 'absolute');
    fprintf('Reflectance error %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(emEst, emRef, 'absolute');
    fprintf('Emission error %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(emEst, emRef, 'normalized');
    fprintf('Emission error (normalized) %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(exEst, exRef, 'normalized');
    fprintf('Excitation error %.3f, std %.3f\n',err,std);
    
end


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



