% Use simulation data to perform tuing parameter crossvalidation for the
% multi fluorophore algorithm (parameters alpha, beta and eta). Note that this
% function is computation heavy and when not parallelized can take a lot of
% time to complete.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Choose save directory. If empty no results will be saved.
dirName = fullfile(fiToolboxRootPath,'results','xVal');
if ~exist(dirName,'dir'), mkdir(dirName); end
% saveFName = fullfile(dirName,[inFName '_xVal_multiFl.mat']);
saveFName = [];

% Load simulation data
% The simulation data file contains wavelength sampling vector.
inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

maxIter = 250;                              % Number of algorithm iterations                   

alphaRange = logspace(-3,1,10);         % Parameter values
betaRange = logspace(-3,1,10);
etaRange = logspace(-3,1,10);

nAlpha = length(alphaRange);
nBeta = length(betaRange);
nEta = length(etaRange);

[alpha, beta, eta] = ndgrid(alphaRange, betaRange, etaRange);
nSets = numel(alpha);

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
       
nSamples = size(measVals,3);
cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);


totalPixelErr = zeros(nAlpha,nBeta,nEta);
reflPixelErr = zeros(nAlpha,nBeta,nEta);
flPixelErr = zeros(nAlpha,nBeta,nEta);
reflErr = zeros(nAlpha,nBeta,nEta);
dMatErr = zeros(nAlpha,nBeta,nEta);
dMatNormErr = zeros(nAlpha,nBeta,nEta);

totalPixelStd = zeros(nAlpha,nBeta,nEta);
reflPixelStd = zeros(nAlpha,nBeta,nEta);
flPixelStd = zeros(nAlpha,nBeta,nEta);
reflStd = zeros(nAlpha,nBeta,nEta);
dMatStd = zeros(nAlpha,nBeta,nEta);
dMatNormStd = zeros(nAlpha,nBeta,nEta);


%% The main cross-validation loop

try
    parpool open local
catch 
end

parfor i=1:nSets
    
    [ reflEst, reflCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMatEst, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndMultiFl( measVals, camera, illuminant, cameraGain*deltaL,...
        cameraOffset, reflBasis, emBasis, exBasis, alpha(i), beta(i), beta(i), eta(i), 'maxIter',maxIter);
    
    % Compute errors
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [totalPixelErr(i), totalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals,[nChannels*nFilters,nSamples]), 'absolute');
    [reflPixelErr(i), reflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), 'absolute');
    [flPixelErr(i), flPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), 'absolute');
    
    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, reflRef, 'absolute');
    
    [dMatErr(i), dMatStd(i)] = fiComputeError(dMatEst, dMatRef, 'absolute');
    [dMatNormErr(i), dMatNormStd(i)] = fiComputeError(dMatEst, dMatRef, 'normalized');
    
end

try
    parpool close
catch
end


%% Save results

if ~isempty(saveFName)
    save(saveFName,'alpha','beta','eta','alphaRange','betaRange','etaRange','nAlpha','nBeta','nEta',...
        'totalPixelErr','reflPixelErr','flPixelErr','reflErr','dMatErr','dMatNormErr',...
        'totalPixelStd','reflPixelStd','flPixelStd','reflStd','dMatStd','dMatNormStd');
end

