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

alphaRange = logspace(-2,1,3);
betaRange = logspace(-2,1,3);
nuRange = logspace(-2,1,3);

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
    [dMatNormErr(i), dMatNorm(i)] = fiComputeError(dMatEst, dMatRef, 'normalized');
    
end

try
    matlabpool close
catch
end

%% Save results
dirName = fullfile(fiToolboxRootPath,'results');
if ~exist(dirName,'dir'), mkdir(dirName); end

fName = fullfile(dirName,[inFName '_xVal.mat']);

save(fName,'alpha','beta','gamma',...
           'totalPixelErr','reflPixelErr','flPixelErr','reflErr','dMatErr','dMatNormErr',...
           'totalPixelStd','reflPixelStd','flPixelStd','reflStd','dMatStd','dMatNormStd');


