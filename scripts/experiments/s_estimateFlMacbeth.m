close all;
clear variables;
clc;

ieInit;

%%

nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

alpha = 0.0;
beta = 0.0;

testFileName = 'Macbeth+Fl';
backgroundFileName = 'Background';

wave = 380:4:1000;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Create basis function sets
reflBasis = createBasisSet('reflectance','wave',wave','n',nReflBasis);
exBasis = createBasisSet('excitation','wave',wave','n',nExBasis);
emBasis = createBasisSet('emission','wave',wave','n',nEmBasis);


% Load the light spectra (in photons). We scale by the maximum for
% numerical stability. This does not matter in the long run as we are
% calibrating for an unknown gain paramter.
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
illuminantPhotons = Energy2Quanta(wave,illuminant);
illuminantPhotons = illuminantPhotons/max(illuminantPhotons(:));

nChannels = size(illuminant,2);

% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);


% Load the calibration target reflectance
% fName = fullfile(fiToolboxRootPath,'data','experiments','chalk');
% calibRefl = ieReadSpectra(fName,wave);
% calibRefl = ones(nWaves,1);

% Load the test target reflectance
fName = fullfile(fiToolboxRootPath,'data','experiments','macbethChart');
reflRef = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'data','redSlideTransmittance');
redTr = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'data','greenSlideTransmittance');
greenTr = ieReadSpectra(fName,wave);

reflRef(:,[1:4:24 2:4:24]) = diag(greenTr)*reflRef(:,[1:4:24 2:4:24]);
reflRef(:,[3:4:24 4:4:24]) = diag(redTr)*reflRef(:,[3:4:24 4:4:24]);


% Load fluorescence data and get reference spectra
fName = fullfile(fiToolboxRootPath,'data','redSlideFluorescence');
redFl = fiReadFluorophore(fName,'wave',wave);

fName = fullfile(fiToolboxRootPath,'data','greenSlideFluorescence');
greenFl = fiReadFluorophore(fName,'wave',wave);

fluorophores = [greenFl; redFl];
flScene = fluorescentSceneCreate('type','fromfluorophore','fluorophore',fluorophores,'wave',wave);

dMatRef = fluorescentSceneGet(flScene,'donaldsonReference','sceneSize',[4 6]);
exRef = fluorescentSceneGet(flScene,'excitationReference','sceneSize',[4 6]);
emRef = fluorescentSceneGet(flScene,'emissionReference','sceneSize',[4 6]);

%% Calibration
% The whole linear image formation model has an unknown gain parameter, we
% are computing this gain for every pixel by taking an image of a white
% surface with a known reflectance. The gain is just the result of the
% division between pixel intensities and linear model predictions.

% Predict the response of a target surface
prediction = deltaL*((camera')*illuminantPhotons);

% Generate the gain map for every pixel
fName = fullfile(fiToolboxRootPath,'data','experiments',backgroundFileName);
[~, ~, scaledRAW] = fiReadImageStack(fName);
hh = size(scaledRAW,1);
ww = size(scaledRAW,2);
gainMap = repmat(shiftdim(prediction,-2),[hh ww 1 1])./scaledRAW;


%% Extract data from a Macbeth image
fName = fullfile(fiToolboxRootPath,'data','experiments',testFileName);
[RAW, ~, scaledMacbeth] = fiReadImageStack(fName);
linearVals = scaledMacbeth.*gainMap;


% Read the sensor data
cp = [35 875;1246 926;1278 134;67 68];
measVals = zeros(nFilters,nChannels,24);
for f=1:nFilters
    sensor = createCameraModel(f);

    for i=1:nChannels
       
        sensor = sensorSet(sensor,'volts',linearVals(:,:,f,i));
        ieAddObject(sensor);

        tmp = macbethSelect(sensor,0,1,cp);
        measVals(f,i,:) = cellfun(@(x) mean(x(~isnan(x) & ~isinf(x))),tmp);
    end 
end

nSamples = size(measVals,3);

%% Estimate the reflectance and single fluorescence

cameraGain = ones([nFilters, nChannels, 24]);
cameraOffset = zeros([nFilters, nChannels, 24]);


[ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
fiRecReflAndFl( measVals, camera, cameraGain*deltaL, cameraOffset, illuminantPhotons, reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter',20);


measValsEst = reflValsEst + flValsEst;


[err, std] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,24]), reshape(measVals - cameraOffset,[nChannels*nFilters,24]), 'default');
fprintf('Total pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,24]), reshape(reflValsRef,[nChannels*nFilters,24]), 'default');
fprintf('Reflected pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,24]), reshape(flValsRef,[nChannels*nFilters,24]), 'default');
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


% Prediction for different filters
tmp1 = reshape(measValsEst,[nFilters nChannels*nSamples])';
tmp2 = reshape(measVals,[nFilters nChannels*nSamples])';

figure;
hold all; grid on; box on;
plot(tmp1(:,:),tmp2(:,:),'.');
xlabel('Model predicted pixel value');
ylabel('ISET pixel value');

% Prediction for different channels
tmp1 = reshape(permute(measValsEst,[2 3 1]),[nChannels nFilters*nSamples]);
tmp2 = reshape(permute(measVals,[2 3 1]),[nChannels nFilters*nSamples]);
figure;
hold all; grid on; box on;
plot(tmp1',tmp2','.');
xlabel('Model predicted pixel value');
ylabel('ISET pixel value');


% Pixel prediction per patch
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    
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
    ylim([0.0 1.05*max(emRef(:))]);

    rmse = sqrt(mean((emEst(:,sampleID) - emRef(:,sampleID)).^2));
    title(sprintf('RMSE %.2f',rmse));

end
end



