close all;
clear all;
clc;

inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

alpha = 0.01;
gamma = 0.1;

% Create basis function sets
nReflBasis = 5;
nExBasis = 12;

[reflBasis, reflScore] = createBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = createBasisSet('excitation','wave',wave','n',nExBasis);


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

% Create a database of emission spectra
fName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell');
[flSet, fluorophoreIDs] = fiReadFluorophoreSet(fName,'wave',wave,...
            'peakEmRange',[wave(5) wave(end-5)],...
            'peakExRange',[wave(5) wave(end-5)]);

tmpScene = fluorescentSceneCreate('type','fromfluorophore','fluorophore',flSet,'wave',wave);
DB = fluorescentSceneGet(tmpScene,'emissionReference');

       
%% Load simulation data 


nSamples = size(measVals,3);
cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

%%

[ reflEst, reflCoeffs, emEst, emChromaticity, exEst, exCoeffs, reflValsEst, flValsEst, hist ] = fiRecReflAndFlMultistep( measVals,...
    camera, cameraGain*deltaL, cameraOffset, illuminant, reflBasis, DB, exBasis, alpha, gamma );

measValsEst = reflValsEst + flValsEst;


[err, std] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals - cameraOffset,[nChannels*nFilters,nSamples]), '');
fprintf('Total pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), '');
fprintf('Reflected pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), '');
fprintf('Fluoresced pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reflEst, reflRef, '');
fprintf('Reflectance error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(emEst, emRef, '');
fprintf('Emission error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(emEst, emRef, 'normalized');
fprintf('Emission error (normalized) %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(exEst, exRef, 'normalized');
fprintf('Excitation error %.3f, std %.3f\n',err,std);



%% Plot the results

% Predicted vs. simulated pixel intensities

figure;
hold all; grid on; box on;
plot(measValsEst(:),measVals(:),'.');
xlabel('Model predicted pixel value');
ylabel('ISET pixel value');

figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on; axis square;

    tmp1 = measValsEst(:,:,sampleID);
    tmp2 = measVals(:,:,sampleID);
    plot(tmp1(:),tmp2(:),'.');

end
end


% Convergence
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot([hist{sampleID}.objValsRefl, hist{sampleID}.objValsChr],'LineWidth',2);
    

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


