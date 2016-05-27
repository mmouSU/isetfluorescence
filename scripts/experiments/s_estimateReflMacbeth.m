% This script uses a reflectance estimation algorithm (fiRecRefl) to
% estimate surface spectral reflectance (no fluorescence) of a Macbeth test
% chart.  
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

ieInit;

%%

nSamples = 24;

nBasis = 12;
testFileName = 'Macbeth';
backgroundFileName = 'Background';

wave = 380:4:1068;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Create basis function sets
nReflBasis = 5;
reflBasis = fiCreateBasisSet('reflectance','wave',wave','n',nReflBasis);


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

% Load the test target reflectance
fName = fullfile(fiToolboxRootPath,'data','macbethChart');
reflRef = ieReadSpectra(fName,wave);


%% Calibration
% The whole linear image formation model has an unknown gain parameter, we
% are computing this gain for every pixel by taking an image of a white
% surface with a known reflectance. The gain is just the result of the
% division between pixel intensities and linear model predictions.

% Predict the response of a target surface
prediction = deltaL*((camera')*illuminantPhotons);

% Generate the gain map for every pixel
fName = fullfile(fiToolboxRootPath,'data','experiments',backgroundFileName);
[~, ~, scaledRAW, shutterBackground] = fiReadImageStack(fName);
hh = size(scaledRAW,1);
ww = size(scaledRAW,2);

% For every illuminant channel and the monochromatic filter pick a
% reference point (say in the middle). Compute how would you need to scale
% other pixels so that the image intensity is uniform. We are using the
% monochromatic channel for the computations and then using the same map
% across all the channels
refPt = scaledRAW(hh/2,ww/2,1,:);
refImg = repmat(refPt,[hh ww 1 1]);
scaleMap = refImg./scaledRAW(:,:,1,:);
scaleMap = repmat(scaleMap,[1 1 nFilters 1]);

nonuniformityCorrected = scaledRAW.*scaleMap;

% Now for each channel compute the gain parameter, that represents the
% precise value of the light intensity.
avgIntensity = squeeze(mean(mean(nonuniformityCorrected(:,:,1,:),1),2));
gains = avgIntensity'./prediction(1,:);

% Assume that the light intensity is constant, we can apply the same gain
% to different camera filters.
cameraGain = repmat(gains,[nFilters, 1, nSamples]);
cameraOffset = zeros([nFilters, nChannels, nSamples]);


%% Extract data from a Macbeth image
fName = fullfile(fiToolboxRootPath,'data','experiments',testFileName);
[~, ~, scaledMacbeth] = fiReadImageStack(fName);
linearVals = scaledMacbeth.*scaleMap;


% Read the sensor data
cp = [35 875;1246 926;1278 134;67 68];
measVals = zeros(nFilters,nChannels,nSamples);
for f=1:nFilters
    sensor = createCameraModel(f);

    for i=1:nChannels
       
        sensor = sensorSet(sensor,'volts',linearVals(:,:,f,i));
        ieAddObject(sensor);

        tmp = macbethSelect(sensor,0,1,cp);
        measVals(f,i,:) = cellfun(@(x) mean(x(~isnan(x) & ~isinf(x))),tmp);
    end 
end

% Normalize the measured pixel intensities, so that the maxium for each 
% patch is 1. To preserve the image formation model we need to scale camera
% gains accordingly.
nF = max(max(measVals,[],1),[],2);
nF = repmat(nF,[nFilters nChannels 1]);
measVals = measVals./nF;
cameraGain = cameraGain./nF;

%% Estimate the reflectance
% and cross-validate for the most optimal smoothness parameter setting.

alphaSet = logspace(-2,0,100);

[optAlpha, rmsError, reflEst, reflCoeffs, predVals] = fiXValRefl(measVals,...
    camera,...
    cameraGain*deltaL,...
    cameraOffset,...
    illuminantPhotons,...
    reflBasis,...
    alphaSet,...
    reflRef);

figure; 
plot(alphaSet,rmsError);
xlabel('\alpha');
ylabel('RMS error');
title(sprintf('Optimal \alpha = %f',optAlpha));


% Measured vs. predicted pixel intensities
figure;
hold on; grid on; box on;
plot(measVals(:) - cameraOffset(:),predVals(:),'.');
xlabel('Measured pixel intensity');
ylabel('Model prediction');


% Estimated reflectance
figure;
for x=1:6
    for y=1:4
        
        indx = y + (x-1)*4;
    
        subfigIndx = (y-1)*6 + x;
        
        subplot(4,6,subfigIndx);
        hold on; box on; grid on;
        plot(wave,reflEst(:,indx),'g','LineWidth',2);
        plot(wave,reflRef(:,indx),'b--','LineWidth',2);

        rmse = sqrt(mean(reflEst(:,indx) - reflRef(:,indx)).^2);
        title(sprintf('RMSE %.2f',rmse));
        
        ylim([0 1]);
        xlim([min(wave), max(wave)]);
        
    end
end


