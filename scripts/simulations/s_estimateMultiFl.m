% Use the multi-fluorophore algorithm to estimate the reflectance and fluorescence
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

% Wavelength sampling is specified in the data file.
deltaL = wave(2) - wave(1);
nWaves = length(wave);

maxIter = 25;

alpha = 0.1;
beta = 0.1;
eta = 0.1;

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
cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

[ reflEst, reflCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMatEst, reflValsEst, flValsEst, hist  ] = ...
    fiRecReflAndMultiFl( measVals, camera, illuminant, cameraGain*deltaL,...
                         cameraOffset, reflBasis, emBasis, exBasis, alpha, beta, beta, eta, 'maxIter',maxIter,...
                         'dMatRef',dMatRef,'reflRef',reflRef,'pixelRef',true);



%% Compute errors

measValsEst = reflValsEst + flValsEst + cameraOffset;

[err, std] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals,[nChannels*nFilters,nSamples]), 'absolute');
fprintf('Total pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), 'absolute');
fprintf('Reflected pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), 'absolute');
fprintf('Fluoresced pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reflEst, reflRef, 'absolute');
fprintf('Reflectance error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(dMatEst, dMatRef, 'absolute');
fprintf('Donaldson Matrix error %f, std %f\n',err,std);

[err, std] = fiComputeError(dMatEst, dMatRef, 'normalized');
fprintf('Donaldson Matrix (normalized) %.3f, std %.3f\n',err,std);




%% Plot the results


figure;
hold all; grid on; box on;
plot(measValsEst(:),measVals(:),'.');
xlim([0 1]);
ylim([0 1]);
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
    plot([hist{sampleID}.prRes, hist{sampleID}.dualRes],'LineWidth',2);
    xlim([0 length(hist{sampleID}.prRes)]);
    ylim([1e-5 10]);
    set(gca,'yscale','log');
end
end

% Estimate convergence
pixelErr = zeros(maxIter,1);
dMatErr = zeros(maxIter,1);
reflErr = zeros(maxIter,1);
for i=1:nSamples
    pixelErr = pixelErr + hist{i}.pixelErr;
    dMatErr = dMatErr + hist{i}.dMatErr;
    reflErr = reflErr + hist{i}.reflErr;
end

figure; 
hold on; grid on; box on;
plot([pixelErr, dMatErr, reflErr]/nSamples);
xlim([1 maxIter]);
set(gca,'yscale','log');
set(gca,'xscale','log');
legend({'Pixel','Donaldson','Reflectance'});


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

% Estimated vs. ground truth Donaldson matrices: scatter plot
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;

    plot(dMatEst{sampleID}(:),dMatRef{sampleID}(:),'.');

    rmse = sqrt(mean((dMatEst{sampleID}(:) - dMatRef{sampleID}(:)).^2));
    title(sprintf('RMSE %.2e',rmse));
end
end

% Estimated vs. ground truth Donaldson matrics: shape
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    
    data = [dMatEst{sampleID} dMatRef{sampleID}];
    imagesc(data);

end
end
