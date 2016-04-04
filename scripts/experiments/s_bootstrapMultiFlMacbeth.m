close all;
clear variables;
clc;

ieInit;

%%

nBootstrap = 100;
sampleSize = 1;

nSamples = 24;

nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;


alpha = 2;
beta = 5;
eta = 0.1;

testFileName = 'Macbeth+multiFl2';
backgroundFileName = 'Background+multiFl';

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

fName = fullfile(fiToolboxRootPath,'data','redFlTransmittance');
redTr = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'data','greenFlTransmittance');
greenTr = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'data','amberFlTransmittance');
amberTr = ieReadSpectra(fName,wave);

reflRef(:,1:4:24) = diag(amberTr)*reflRef(:,1:4:24);
reflRef(:,2:4:24) = diag(greenTr)*diag(amberTr)*reflRef(:,2:4:24);
reflRef(:,3:4:24) = diag(greenTr)*diag(redTr)*reflRef(:,3:4:24);
reflRef(:,4:4:24) = diag(redTr)*reflRef(:,4:4:24);

% Load fluorescence data and get reference spectra
fName = fullfile(fiToolboxRootPath,'data','redFl');
redFl = fiReadFluorophore(fName,'wave',wave);

fName = fullfile(fiToolboxRootPath,'data','greenFl');
greenFl = fiReadFluorophore(fName,'wave',wave);

fName = fullfile(fiToolboxRootPath,'data','amberFl');
amberFl = fiReadFluorophore(fName,'wave',wave);


fluorophores(1,1,1) = amberFl;
fluorophores(1,1,2) = amberFl;

fluorophores(2,1,1) = greenFl;
fluorophores(2,1,2) = amberFl;

fluorophores(3,1,1) = greenFl;
fluorophores(3,1,2) = redFl;

fluorophores(4,1,1) = redFl;
fluorophores(4,1,2) = redFl;


flQe(1,1,1) = 0.25/2;
flQe(1,1,2) = 0.25/2;

flQe(2,1,1) = 0.20;
flQe(2,1,2) = 0.14;

flQe(3,1,1) = 0.20;
flQe(3,1,2) = 0.14;

flQe(4,1,1) = 0.23/2;
flQe(4,1,2) = 0.23/2;

% fluorophores = [greenFl; redFl];
flScene = fluorescentSceneCreate('type','fromfluorophore','fluorophore',fluorophores,'wave',wave,'qe',flQe);

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
[RAW, ~, scaledRAW, shutterBackground] = fiReadImageStack(fName);
hh = size(scaledRAW,1);
ww = size(scaledRAW,2);

%{
for f=1:nFilters
    figure;
    for c=1:nChannels
        subplot(4,4,c);
        imagesc(RAW(:,:,f,c));
    end
end
%}

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
[reflValsRef, flValsRef] = fiComputeReflFlContrib(camera,illuminantPhotons,cameraGain(:,:,1)*deltaL,reflRef,dMatRef);


%% Extract data from a Macbeth image
fName = fullfile(fiToolboxRootPath,'data','experiments',testFileName);
[RAW, ~, scaledMacbeth, shutterMacbeth] = fiReadImageStack(fName);
linearVals = scaledMacbeth.*scaleMap;

% Read the sensor data
cp = [51 846;1275 873;1278 81;58 61];
for f=1:nFilters
    sensor = createCameraModel(f);

    for i=1:nChannels
       
        sensor = sensorSet(sensor,'volts',linearVals(:,:,f,i));
        % sensor = sensorSet(sensor,'volts',RAW(:,:,f,i));
        ieAddObject(sensor);

        [tmp, ~, ~, cp] = macbethSelect(sensor,0,1,cp);

        if f==1 && i==1
            nPixels = length(tmp{1});
            measVals = zeros(nFilters,nChannels,24,nPixels);
        end
        
        measVals(f,i,:,:) = cell2mat(tmp)';
    end 
end


% Normalize the measured pixel intensities, so that the maxium for each 
% patch is 1. To preserve the image formation model we need to scale camera
% gains accordingly.
nF = max(max(mean(measVals,4),[],1),[],2);
nF = repmat(nF,[nFilters nChannels 1]);
measVals = measVals./repmat(nF,[1 1 1 nPixels]);
cameraGain = cameraGain./nF;


%% Estimate the reflectance and single fluorescence


reflEst = cell(nBootstrap,1);
rfCoeffs = cell(nBootstrap,1);
emEst = cell(nBootstrap,1);
emCoeffs = cell(nBootstrap,1);
exEst = cell(nBootstrap,1);
exCoeffs = cell(nBootstrap,1);
dMatEst = cell(nBootstrap,1);
reflValsEst = cell(nBootstrap,1);
flValsEst = cell(nBootstrap,1);
hist = cell(nBootstrap,1);
inds = cell(nBootstrap,1);

try 
    matlabpool open local
catch
end

parfor b=1:nBootstrap
    
    % Pick a bootstrap sample, and average selected pixels
    inds{b} = randi(nPixels,[sampleSize 1]);
    bootstrapMeasVals = mean(measVals(:,:,:,inds{b}),4);
    

    [ reflEst{b}, reflCoeffs{b}, emEst{b}, emCoeffs{b}, exEst{b}, exCoeffs{b}, dMatEst{b}, reflValsEst{b}, flValsEst{b}, hist{b}  ] = ...
    fiRecReflAndMultiFl( bootstrapMeasVals, camera, illuminantPhotons, cameraGain*deltaL,...
                         cameraOffset, reflBasis, emBasis, exBasis, alpha, beta, beta, eta, 'maxIter',500,'rescaleRho',false);

end

try
    matlabbool close
catch 
end


dirName = fullfile(fiToolboxRootPath,'results','bootstrap');
if ~exist(dirName,'dir'), mkdir(dirName); end;

fName = fullfile(dirName,sprintf('multifl_%s_bootstrap_%i_%i.mat',testFileName,nBootstrap,sampleSize));
save(fName,'reflEst','reflCoeffs','emEst','emCoeffs','exEst','exCoeffs','dMatEst','reflValsEst','flValsEst','hist',...
            'wave','alpha','beta','eta','reflRef','exRef','emRef','dMatRef','measVals','nBootstrap','sampleSize');
                     
                     
                     
