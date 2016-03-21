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
nNu = 10;

alphaRange = logspace(-3,1,nAlpha);
betaRange = logspace(-3,1,nBeta);
nuRange = logspace(-3,1,nNu);

[alpha, beta, nu] = ndgrid(alphaRange, betaRange, nuRange);
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
dMatErr = zeros(nSets,1);
dMatNormErr = zeros(nSets,1);

totalPixelStd = zeros(nSets,1);
reflPixelStd = zeros(nSets,1);
flPixelStd = zeros(nSets,1);
reflStd = zeros(nSets,1);
dMatStd = zeros(nSets,1);
dMatNormStd = zeros(nSets,1);


%% The main cross-validation loop

try
    matlabpool open local
catch 
end

parfor i=1:nSets
    
    [ reflEst, reflCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMatEst, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndMultiFl( measVals, camera, illuminant, cameraGain*deltaL,...
        cameraOffset, reflBasis, emBasis, exBasis, alpha(i), beta(i), beta(i), nu(i), 'maxIter',250);
    
    
    % Compute errors
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [totalPixelErr(i), totalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals,[nChannels*nFilters,nSamples]), '');
    [reflPixelErr(i), reflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), '');
    [flPixelErr(i), flPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), '');
    
    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, reflRef, '');
    
    [dMatErr(i), dMatStd(i)] = fiComputeError(dMatEst, dMatRef, '');
    [dMatNormErr(i), dMatNormStd(i)] = fiComputeError(dMatEst, dMatRef, 'normalized');
    
end

try
    matlabpool close
catch
end

totalPixelErr = reshape(totalPixelErr,[nAlpha, nBeta, nNu]);
reflPixelErr = reshape(reflPixelErr,[nAlpha, nBeta, nNu]);
flPixelErr = reshape(flPixelErr,[nAlpha, nBeta, nNu]);
reflErr = reshape(reflErr,[nAlpha, nBeta, nNu]);
dMatErr = reshape(dMatErr,[nAlpha, nBeta, nNu]);
dMatNormErr = reshape(dMatNormErr,[nAlpha, nBeta, nNu]);

totalPixelStd = reshape(totalPixelStd,[nAlpha, nBeta, nNu]);
reflPixelStd = reshape(reflPixelStd,[nAlpha, nBeta, nNu]);
flPixelStd = reshape(flPixelStd,[nAlpha, nBeta, nNu]);
reflStd = reshape(reflStd,[nAlpha, nBeta, nNu]);
dMatStd = reshape(dMatStd,[nAlpha, nBeta, nNu]);
dMatNormStd = reshape(dMatNormStd,[nAlpha, nBeta, nNu]);



%% Save results

dirName = fullfile(fiToolboxRootPath,'results','xVal');
if ~exist(dirName,'dir'), mkdir(dirName); end

fName = fullfile(dirName,[inFName '_xVal_multiFl.mat']);

save(fName,'alpha','beta','nu','alphaRange','betaRange','nuRange',...
           'totalPixelErr','reflPixelErr','flPixelErr','reflErr','dMatErr','dMatNormErr',...
           'totalPixelStd','reflPixelStd','flPixelStd','reflStd','dMatStd','dMatNormStd');


%% Display results

lineStyle = {'rs-','gd-','bo-'};
lw = 2;

% Alpha

[pixelErrPlot, minLoc] = min(reshape(permute(totalPixelErr,[2 3 1]),[nBeta*nNu, nAlpha]));
tmp = reshape(permute(totalPixelStd,[2 3 1]),[nBeta*nNu, nAlpha]);
inds = sub2ind([nBeta*nNu, nAlpha],minLoc,1:nAlpha);
pixelErrPlotStd = tmp(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reshape(permute(reflErr,[2 3 1]),[nBeta*nNu, nAlpha]));
tmp = reshape(permute(reflStd,[2 3 1]),[nBeta*nNu, nAlpha]);
inds = sub2ind([nBeta*nNu, nAlpha],minLoc,1:nAlpha);
reflErrPlotStd = tmp(inds)/sqrt(24);

[dMatNormErrPlot, minLoc] = min(reshape(permute(dMatNormErr,[2 3 1]),[nBeta*nNu, nAlpha]));
tmp = reshape(permute(dMatNormStd,[2 3 1]),[nBeta*nNu, nAlpha]);
inds = sub2ind([nBeta*nNu, nAlpha],minLoc,1:nAlpha);
dMatNormPlotStd = tmp(inds)/sqrt(24);

figure;
hold all; grid on; box on;
pl(1) = errorbar(alphaRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(alphaRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(alphaRange,dMatNormErrPlot,dMatNormPlotStd,lineStyle{3});
set(pl,'lineWidth',lw);
set(gca,'xscale','log');
ylim([0 0.05]);
xlim([0.95*min(alphaRange) 1.05*max(alphaRange)]);
xlabel('\alpha');
ylabel('min_{\beta,\nu}(RMSE)');

% Beta

[pixelErrPlot, minLoc] = min(reshape(permute(totalPixelErr,[1 3 2]),[nAlpha*nNu, nBeta]));
tmp = reshape(permute(totalPixelStd,[1 3 2]),[nAlpha*nNu, nBeta]);
inds = sub2ind([nAlpha*nNu, nBeta],minLoc,1:nBeta);
pixelErrPlotStd = tmp(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reshape(permute(reflErr,[1 3 2]),[nAlpha*nNu, nBeta]));
tmp = reshape(permute(reflStd,[1 3 2]),[nAlpha*nNu, nBeta]);
inds = sub2ind([nAlpha*nNu, nBeta],minLoc,1:nBeta);
reflErrPlotStd = tmp(inds)/sqrt(24);

[dMatNormErrPlot, minLoc] = min(reshape(permute(dMatNormErr,[1 3 2]),[nAlpha*nNu, nBeta]));
tmp = reshape(permute(dMatNormStd,[1 3 2]),[nAlpha*nNu, nBeta]);
inds = sub2ind([nAlpha*nNu, nBeta],minLoc,1:nBeta);
dMatNormPlotStd = tmp(inds)/sqrt(24);



figure;
hold all; grid on; box on;
pl(1) = errorbar(betaRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(betaRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(betaRange,dMatNormErrPlot,dMatNormPlotStd,lineStyle{3});
set(gca,'xscale','log');
set(pl,'lineWidth',lw);
xlim([0.95*min(betaRange) 1.05*max(betaRange)]);
ylim([0 0.05]);
xlabel('\beta');
ylabel('min_{\alpha,\nu}(RMSE)');


% Nu

[pixelErrPlot, minLoc] = min(reshape(permute(totalPixelErr,[1 2 3]),[nAlpha*nBeta, nNu]));
tmp = reshape(permute(totalPixelStd,[1 2 3]),[nAlpha*nBeta, nNu]);
inds = sub2ind([nAlpha*nBeta, nNu],minLoc,1:nNu);
pixelErrPlotStd = tmp(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reshape(permute(reflErr,[1 2 3]),[nAlpha*nBeta, nNu]));
tmp = reshape(permute(reflStd,[1 2 3]),[nAlpha*nBeta, nNu]);
inds = sub2ind([nAlpha*nBeta, nNu],minLoc,1:nNu);
reflErrPlotStd = tmp(inds)/sqrt(24);

[dMatNormErrPlot, minLoc] = min(reshape(permute(dMatNormErr,[1 2 3]),[nAlpha*nBeta, nNu]));
tmp = reshape(permute(dMatNormStd,[1 2 3]),[nAlpha*nBeta, nNu]);
inds = sub2ind([nAlpha*nBeta, nNu],minLoc,1:nNu);
dMatNormPlotStd = tmp(inds)/sqrt(24);



figure;
hold all; grid on; box on;
pl(1) = errorbar(nuRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(nuRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(nuRange,dMatNormErrPlot,dMatNormPlotStd,lineStyle{3});
set(gca,'xscale','log');
set(pl,'lineWidth',lw);
xlim([0.95*min(nuRange) 1.05*max(nuRange)]);
ylim([0 0.05]);
xlabel('\nu');
ylabel('min_{\alpha,\beta}(RMSE)');
