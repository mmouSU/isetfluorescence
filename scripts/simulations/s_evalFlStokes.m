% Evaluate the single fluorophore algorithm on fluorophores with different
% values of the Stokes shift.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Save results to a file saveFName if not empty.
dirName = fullfile(fiToolboxRootPath,'results','evaluation');
if ~exist(dirName,'dir'), mkdir(dirName); end
% saveFName = fullfile(dirName,[dataset '_simStokes_Fl.mat']);
saveFName = [];

wave = 380:4:1000;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Tuning parameters
alpha = 0.1;
beta = 0.1;

% Scene properties
dataset = 'McNamara-Boswell';
height = 4;
width = 6;
nSamples = height*width;
nFluorophores = 1;
flQe = 0.1;


% Stokes shift ranges
stokesRanges = [0 25 50 75 100];
nStokes = length(stokesRanges)-1;

% Create basis function sets
nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = fiCreateBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = fiCreateBasisSet('excitation','wave',wave','n',nExBasis);
[emBasis, emScore] = fiCreateBasisSet('emission','wave',wave','n',nEmBasis);

% Error placeholder variables
totalPixelErr = zeros(nStokes,1);
reflPixelErr = zeros(nStokes,1);
flPixelErr = zeros(nStokes,1);
reflErr = zeros(nStokes,1);
exErr = zeros(nStokes,1);
exNormErr = zeros(nStokes,1);
emErr = zeros(nStokes,1);
emNormErr = zeros(nStokes,1);

totalPixelStd = zeros(nStokes,1);
reflPixelStd = zeros(nStokes,1);
flPixelStd = zeros(nStokes,1);
reflStd = zeros(nStokes,1);
exStd = zeros(nStokes,1);
exNormStd = zeros(nStokes,1);
emStd = zeros(nStokes,1);
emNormStd = zeros(nStokes,1);

%% The main cross-validation loop

try
    cluster = parcluster('local');
    cluster.NumWorkers = min(nStokes,35);
    pool = parpool(cluster,cluster.NumWorkers);
catch 
end

parfor i=1:nStokes

    fName = fullfile(fiToolboxRootPath,'data','simulations',sprintf('%s_%ix%ix%i_stokes_%i.mat',dataset,height,width,nFluorophores,stokesRanges(i+1)));
    data = load(fName);

    
    cameraGain = repmat(data.cameraGain,[1 1 nSamples]);
    cameraOffset = repmat(data.cameraOffset,[1 1 nSamples]);
    
    
    [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndFl( data.measVals, data.camera, cameraGain*deltaL, cameraOffset, data.illuminantPhotons,...
        reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter', 250 );

    
    %% Evaluation
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [totalPixelErr(i), totalPixelStd(i)] = fiComputeError(reshape(measValsEst,[data.nChannels*data.nFilters,nSamples]), reshape(data.measVals,[data.nChannels*data.nFilters,nSamples]), 'absolute');
    [reflPixelErr(i), reflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[data.nChannels*data.nFilters,nSamples]), reshape(data.reflValsRef,[data.nChannels*data.nFilters,nSamples]), 'absolute');
    [flPixelErr(i), flPixelStd(i)] = fiComputeError(reshape(flValsEst,[data.nChannels*data.nFilters,nSamples]), reshape(data.flValsRef,[data.nChannels*data.nFilters,nSamples]), 'absolute');
    
    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, data.reflRef, 'absolute');
    
    [emErr(i), emStd(i)] = fiComputeError(emEst, data.emRef, 'absolute');
    [emNormErr(i), emNormStd(i)] = fiComputeError(emEst, data.emRef, 'normalized');
    
    [exErr(i), exStd(i)] = fiComputeError(exEst, data.exRef, 'absolute');
    [exNormErr(i), exNormStd(i)] = fiComputeError(exEst, data.exRef, 'normalized');
    

end

try
    delete(pool);
catch
end

%% Save results

if ~isempty(saveFName)
    save(saveFName,'stokesRanges',...
        'totalPixelErr','reflPixelErr','flPixelErr','reflErr','exErr','exNormErr','emErr','emNormErr',...
        'totalPixelStd','reflPixelStd','flPixelStd','reflStd','exStd','exNormStd','emStd','emNormStd');
end
       
