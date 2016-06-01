% Use simulation data to perform tuing parameter crossvalidation for the
% single fluorophore algorithm (parameters alpha and beta). Note that this
% function is computation heavy and when not parallelized can take a lot of
% time to complete.
%
% Copyright, Henryk Blasinski 2016


close all;
clear all;
clc;

% Define the save file name. If empty results won't be saved
dirName = fullfile(fiToolboxRootPath,'results','xVal');
if ~exist(dirName,'dir'), mkdir(dirName); end
% saveFName = fullfile(dirName,[inFName '_xVal_Fl.mat']);
saveFName = [];

% Load simulation data. The simulation data file contains wavelength sampling vector.
inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

maxIter = 250;                           % Number of single fluorophore algorithm iterations

alphaRange = logspace(-3,1,10);         % Specify the range of parameters
betaRange = logspace(-3,1,10);

nAlpha = length(alphaRange);
nBeta = length(betaRange);

[alpha, beta] = ndgrid(alphaRange, betaRange);
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


% Result placeholder variables
totalPixelErr = zeros(nAlpha,nBeta);
reflPixelErr = zeros(nAlpha,nBeta);
flPixelErr = zeros(nAlpha,nBeta);
reflErr = zeros(nAlpha,nBeta);
exErr = zeros(nAlpha,nBeta);
exNormErr = zeros(nAlpha,nBeta);
emErr = zeros(nAlpha,nBeta);
emNormErr = zeros(nAlpha,nBeta);

totalPixelStd = zeros(nAlpha,nBeta);
reflPixelStd = zeros(nAlpha,nBeta);
flPixelStd = zeros(nAlpha,nBeta);
reflStd = zeros(nAlpha,nBeta);
exStd = zeros(nAlpha,nBeta);
exNormStd = zeros(nAlpha,nBeta);
emStd = zeros(nAlpha,nBeta);
emNormStd = zeros(nAlpha,nBeta);


%% The main cross-validation loop

try
    parpool open local
catch 
end

parfor i=1:nSets
    
    [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
    fiRecReflAndFl( measVals, camera, cameraGain*deltaL, cameraOffset, illuminant,...
            reflBasis, emBasis, exBasis, alpha(i), beta(i), beta(i), 'maxIter', maxIter );
    
    
    % Compute errors
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [totalPixelErr(i), totalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals,[nChannels*nFilters,nSamples]), 'absolute');
    [reflPixelErr(i), reflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), 'absolute');
    [flPixelErr(i), flPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), 'absolute');
    
    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, reflRef, 'absolute');
    
    [emErr(i), emStd(i)] = fiComputeError(emEst, emRef, 'absolute');
    [emNormErr(i), emNormStd(i)] = fiComputeError(emEst, emRef, 'normalized');

    [exErr(i), exStd(i)] = fiComputeError(exEst, exRef, 'absolute');
    [exNormErr(i), exNormStd(i)] = fiComputeError(exEst, exRef, 'normalized');
    
end

try
    parpool close
catch
end

%% Save results

if ~isempty(saveFName)
    save(saveFName,'alpha','beta','alphaRange','betaRange',...
        'totalPixelErr','reflPixelErr','flPixelErr','reflErr','exErr','exNormErr','emErr','emNormErr',...
        'totalPixelStd','reflPixelStd','flPixelStd','reflStd','exStd','exNormStd','emStd','emNormStd');
end

