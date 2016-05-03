% close all;
clear variables;
clc;

ieInit;

%%

nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

dsFac = 10;

alpha = 0.1;
beta = 5;
eta = 0.01;

testFileName = 'Natural';
backgroundFileName = 'Background+Natural';
dirName = fullfile(fiToolboxRootPath,'results','experiments');
if ~exist(dirName,'dir'), mkdir(dirName); end;

wave = 380:4:1000;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Create basis function sets
reflBasis = createBasisSet('reflectance','wave',wave','n',nReflBasis);
exBasis = createBasisSet('excitation','wave',wave','n',nExBasis);
emBasis = createBasisSet('emission','wave',wave','n',nEmBasis);

% exBasis = eye(nWaves);
% emBasis = eye(nWaves);

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


%% Extract data from a Macbeth image
fName = fullfile(fiToolboxRootPath,'data','experiments',testFileName);
[RAW, ~, scaledImg, shutterImg] = fiReadImageStack(fName);
linearVals = scaledImg.*scaleMap;


% Downsample the image
dsLinearVals = linearVals(1:dsFac:end,1:dsFac:end,:,:);


rows = [40];
cols = [8 35 94];

nRows = length(rows);
nCols = length(cols);
nSamples = nRows*nCols;

figure;
imshow(dsLinearVals(:,:,1,1),[],'Border','tight','InitialMagnification',500);
for xx=1:nCols
for yy=1:nRows
    rectangle('Position',[cols(xx) rows(yy) 1 1] ,'edgecolor','red');
end
end

% Assume that the light intensity is constant, we can apply the same gain
% to different camera filters.
cameraGain = repmat(gains,[nFilters, 1, nSamples]);
cameraOffset = zeros([nFilters, nChannels, nSamples]);


    
measVals = reshape(dsLinearVals(rows,cols,:,:),[nSamples, nFilters, nChannels]);
measVals = permute(measVals,[2 3 1]);
    
% Normalize the measured pixel intensities, so that the maxium for each
% patch is 1. To preserve the image formation model we need to scale camera
% gains accordingly.
nF = max(measVals,[],2);
nF = repmat(nF,[1 nChannels 1]);
measVals = measVals./nF;
cameraGainCol = cameraGain./nF;

[ reflEst, ~, emEst, ~, exEst, ~, dMatEst, reflValsEst, flValsEst, hist  ] = ...
    fiRecReflAndMultiFl( measVals(:,2:end,:), camera, illuminantPhotons(:,2:end), cameraGainCol(:,2:end,:)*deltaL,...
    cameraOffset(:,2:end,:), reflBasis, emBasis, exBasis, alpha, beta, beta, eta, 'maxIter',5000,'rescaleRho',false);

%% Results
    
measValsEst = reflValsEst + flValsEst;

fName = fullfile(fiToolboxRootPath,'data','flCmap');
load(fName);

figure;
hold all; grid on; box on;
plot(measValsEst(:),measVals(:),'.');
xlabel('Model predicted pixel value');
ylabel('ISET pixel value');


% Prediction for different filters
tmp1 = reshape(measValsEst,[nFilters nChannels*nSamples])';
tmp2 = reshape(measVals,[nFilters nChannels*nSamples])';
maxVal = max([tmp1(:); tmp2(:)]);

figure;
for f=1:nFilters
    subplot(3,3,f);
    hold all; grid on; box on;
    plot(tmp1(:,f),tmp2(:,f),'.');
    plot(linspace(0,maxVal,10),linspace(0,maxVal,10),'r');
    xlabel('Model');
    ylabel('Measured');
    title(sprintf('F: %i',f));
    xlim([0 maxVal]);
    ylim([0 maxVal]);
end

% Prediction for different chanels
tmp1 = reshape(permute(measValsEst,[2 3 1]),[nChannels nFilters*nSamples])';
tmp2 = reshape(permute(measVals,[2 3 1]),[nChannels nFilters*nSamples])';
maxVal = max([tmp1(:); tmp2(:)]);

figure;
for c=1:nChannels
    subplot(4,4,c);
    hold all; grid on; box on;
    plot(tmp1(:,c),tmp2(:,c),'.');
    plot(linspace(0,maxVal,10),linspace(0,maxVal,10),'r');
    xlabel('Model');
    ylabel('Measured');
    title(sprintf('C: %i',c));
    xlim([0 maxVal]);
    ylim([0 maxVal]);
end

% Convergence
figure;
for xx=1:nCols
for yy=1:nRows

    plotID = (yy-1)*nCols + xx;
    sampleID = (xx-1)*nRows + yy;

    subplot(nRows,nCols,plotID);
    hold all; grid on; box on;
    plot([hist{sampleID}.prRes, hist{sampleID}.dualRes],'LineWidth',2);
    xlim([0 length(hist{sampleID}.prRes)]);
    ylim([1e-5 10]);
    set(gca,'yscale','log');
end
end

% Reflectance
figure;
for xx=1:nCols
for yy=1:nRows

    plotID = (yy-1)*nCols + xx;
    sampleID = (xx-1)*nRows + yy;

    subplot(nRows,nCols,plotID);
    hold all; grid on; box on;
    plot(wave,reflEst(:,sampleID),'g','LineWidth',2);
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);

end
end

% Donaldson matrices:
figure;
set(gcf,'Colormap',flCmap);
for xx=1:nCols
for yy=1:nRows

    plotID = (yy-1)*nCols + xx;
    sampleID = (xx-1)*nRows + yy;

    subplot(nRows,nCols,plotID);
    
    data = [dMatEst{sampleID}];
    imagesc(wave, wave, data);
    axis square;

end
end

