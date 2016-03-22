close all;
clear all;
clc;

%% Load simulation data
% The simulation data file contains wavelength sampling vector.
inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

nAlpha = 10;
nBeta = 10;

alphaRange = logspace(-3,1,nAlpha);
betaRange = logspace(-3,1,nBeta);

[alpha, beta] = ndgrid(alphaRange, betaRange);
nSets = numel(alpha);

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
       

nSamples = size(measVals,3);
cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);



totalPixelErr = zeros(nSets,1);
reflPixelErr = zeros(nSets,1);
flPixelErr = zeros(nSets,1);

reflErr = zeros(nSets,1);

exErr = zeros(nSets,1);
exNormErr = zeros(nSets,1);

emErr = zeros(nSets,1);
emNormErr = zeros(nSets,1);



totalPixelStd = zeros(nSets,1);
reflPixelStd = zeros(nSets,1);
flPixelStd = zeros(nSets,1);

reflStd = zeros(nSets,1);

exStd = zeros(nSets,1);
exNormStd = zeros(nSets,1);

emStd = zeros(nSets,1);
emNormStd = zeros(nSets,1);


%% The main cross-validation loop

try
    matlabpool open local
catch 
end

parfor i=1:nSets
    
    [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
    fiRecReflAndFl( measVals, camera, cameraGain*deltaL, cameraOffset, illuminant,...
            reflBasis, emBasis, exBasis, alpha(i), beta(i), beta(i), 'maxIter', 250 );
    
    
    % Compute errors
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [totalPixelErr(i), totalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals,[nChannels*nFilters,nSamples]), '');
    [reflPixelErr(i), reflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), '');
    [flPixelErr(i), flPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), '');
    
    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, reflRef, '');
    
    [emErr(i), emStd(i)] = fiComputeError(emEst, emRef, '');
    [emNormErr(i), emNormStd(i)] = fiComputeError(emEst, emRef, 'normalized');

    [exErr(i), exStd(i)] = fiComputeError(exEst, exRef, '');
    [exNormErr(i), exNormStd(i)] = fiComputeError(exEst, exRef, 'normalized');
    
end

try
    matlabpool close
catch
end

totalPixelErr = reshape(totalPixelErr,[nAlpha, nBeta]);
reflPixelErr = reshape(reflPixelErr,[nAlpha, nBeta]);
flPixelErr = reshape(flPixelErr,[nAlpha, nBeta]);

reflErr = reshape(reflErr,[nAlpha, nBeta]);

exErr = reshape(exErr,[nAlpha, nBeta]);
exNormErr = reshape(exNormErr,[nAlpha, nBeta]);

emErr = reshape(emErr,[nAlpha, nBeta]);
emNormErr = reshape(emNormErr,[nAlpha, nBeta]);


totalPixelStd = reshape(totalPixelStd,[nAlpha, nBeta]);
reflPixelStd = reshape(reflPixelStd,[nAlpha, nBeta]);
flPixelStd = reshape(flPixelStd,[nAlpha, nBeta]);

reflStd = reshape(reflStd,[nAlpha, nBeta]);

exStd = reshape(exStd,[nAlpha, nBeta]);
exNormStd = reshape(exNormStd,[nAlpha, nBeta]);

emStd = reshape(emStd,[nAlpha, nBeta]);
emNormStd = reshape(emNormStd,[nAlpha, nBeta]);



%% Save results

dirName = fullfile(fiToolboxRootPath,'results','xVal');
if ~exist(dirName,'dir'), mkdir(dirName); end

fName = fullfile(dirName,[inFName '_xVal_Fl.mat']);

save(fName,'alpha','beta','alphaRange','betaRange',...
           'totalPixelErr','reflPixelErr','flPixelErr','reflErr','exErr','exNormErr','emErr','emNormErr',...
           'totalPixelStd','reflPixelStd','flPixelStd','reflStd','exStd','exNormStd','emStd','emNormStd');


%% Display results

lineStyle = {'rs-','gd-','bo-','c^-'};
lw = 2;

% Alpha

[pixelErrPlot, minLoc] = min(totalPixelErr,[],2);
inds = sub2ind([nAlpha, nBeta],1:nAlpha,minLoc);
pixelErrPlotStd = totalPixelErr(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reflErr,[],2);
inds = sub2ind([nAlpha, nBeta],1:nAlpha,minLoc);
reflErrPlotStd = reflErr(inds)/sqrt(24);

[exNormErrPlot, minLoc] = min(exNormErr,[],2);
inds = sub2ind([nAlpha, nBeta],1:nAlpha,minLoc);
exNormPlotStd = exNormErr(inds)/sqrt(24);

[emNormErrPlot, minLoc] = min(emNormErr,[],2);
inds = sub2ind([nAlpha, nBeta],1:nAlpha,minLoc);
emNormPlotStd = emNormErr(inds)/sqrt(24);

figure;
hold all; grid on; box on;
pl(1) = errorbar(alphaRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(alphaRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(alphaRange,exNormErrPlot,exNormPlotStd,lineStyle{3});
pl(4) = errorbar(alphaRange,emNormErrPlot,emNormPlotStd,lineStyle{4});
set(pl,'lineWidth',lw);
set(gca,'xscale','log');
ylim([0 0.5]);
xlim([0.95*min(alphaRange) 1.05*max(alphaRange)]);
xlabel('\alpha');
ylabel('min_{\beta}(RMSE)');

% Beta

[pixelErrPlot, minLoc] = min(totalPixelErr,[],1);
inds = sub2ind([nAlpha, nBeta],minLoc,1:nBeta);
pixelErrPlotStd = totalPixelErr(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reflErr,[],1);
inds = sub2ind([nAlpha, nBeta],minLoc,1:nBeta);
reflErrPlotStd = reflErr(inds)/sqrt(24);

[exNormErrPlot, minLoc] = min(exNormErr,[],1);
inds = sub2ind([nAlpha, nBeta],minLoc,1:nBeta);
exNormPlotStd = exNormErr(inds)/sqrt(24);

[emNormErrPlot, minLoc] = min(emNormErr,[],1);
inds = sub2ind([nAlpha, nBeta],minLoc,1:nBeta);
emNormPlotStd = emNormErr(inds)/sqrt(24);



figure;
hold all; grid on; box on;
pl(1) = errorbar(betaRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(betaRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(betaRange,exNormErrPlot,exNormPlotStd,lineStyle{3});
pl(4) = errorbar(betaRange,emNormErrPlot,emNormPlotStd,lineStyle{4});
set(gca,'xscale','log');
set(pl,'lineWidth',lw);
xlim([0.95*min(betaRange) 1.05*max(betaRange)]);
ylim([0 0.5]);
xlabel('\beta');
ylabel('min_{\alpha}(RMSE)');


