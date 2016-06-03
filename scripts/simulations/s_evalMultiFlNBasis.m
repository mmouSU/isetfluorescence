% Evaluate the multi fluorophore algorithm using simulated data captured
% with a bispectral system to evaluate the number of excitation and
% emission basis used for spectra approximation.
%
% This script is computationally intense and may take a long time to run.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Save results to file if saveFName ~= []
% saveFName = fullfile(fiToolboxRootPath,'results','evaluation','multiFl_nBasis.mat');
saveFName = [];

wave = 380:4:1000;
nWaves = length(wave);
deltaL = wave(2) - wave(1);

% Number of excitation and emission basis to evaluate
maxExBasis = 25;
maxEmBasis = 25;

% Tuning parameters
alpha = 0.001;
beta = 0.001;
eta = 0.001;

% Scene properties
flQe = 0.5;
nFluorophores = 1;
dataset = 'McNamara-Boswell';
height = 4;
width = 6;
nSamples = height*width;

% Reflectance basis
nReflBasis = 5;
[reflBasis, reflScore] = fiCreateBasisSet('reflectance','wave',wave','n',nReflBasis);


% Load the light spectra (in photons)
illuminant = eye(nWaves);
nChannels = size(illuminant,2);

% Load camera spectral properties
camera = eye(nWaves);
nFilters = size(camera,2);

cameraGain = ones(nFilters,nWaves);
cameraOffset = zeros(nFilters,nWaves);

% Create reflective scene
scene = sceneCreate('macbethEE_IR','',wave);
scene = sceneSet(scene,'fov',5);
scene = sceneSet(scene,'distance',1);

fName = fullfile(isetRootPath,'data','surfaces','macbethChart');
reflRef = ieReadSpectra(fName,wave);

% Create fluorescent scene
flScene = fluorescentSceneCreate('height',height,'width',width,'wave',wave,'qe',flQe,'nFluorophores',nFluorophores,...
                                 'peakEmRange',[wave(5) wave(end-5)],...
                                 'peakExRange',[wave(5) wave(end-5)],...
                                 'dataSet',dataset);
         
dMatRef = fluorescentSceneGet(flScene,'Donaldson reference');        


[reflValsRef, flValsRef] = fiComputeReflFlContrib(camera,illuminant,cameraGain*deltaL,reflRef,dMatRef);
measVals = reflValsRef + flValsRef;

cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

nF = max(max(measVals,[],1),[],2);
nF = repmat(nF,[nFilters,nChannels,1]);
measVals = measVals./nF;
cameraGain = cameraGain./nF;

[emBasisGrid, exBasisGrid] = meshgrid(1:maxEmBasis,1:maxExBasis);


% Error placeholder variables
reflErr = zeros(maxEmBasis,maxExBasis);
reflStd = zeros(maxEmBasis,maxExBasis);

dMatErr = zeros(maxEmBasis,maxExBasis);
dMatStd = zeros(maxEmBasis,maxExBasis);

pixelErr = zeros(maxEmBasis,maxExBasis);
pixelStd = zeros(maxEmBasis,maxExBasis);


try
    cluster = parcluster('local');
    cluster.NumWorkers = min(numel(emBasisGrid),35);
    pool = parpool(cluster,cluster.NumWorkers);
catch
end

parfor i=1:numel(emBasisGrid);


    exBasis = createBasisSet('excitation','wave',wave','n',exBasisGrid(i));
    emBasis = createBasisSet('emission','wave',wave','n',emBasisGrid(i));
    
    [ reflEst, ~, emEst, ~, exEst, ~, dMatEst, reflValsEst, flValsEst, hist  ] = ...
    fiRecReflAndMultiFl( measVals, camera, illuminant, cameraGain*deltaL,...
                         cameraOffset, reflBasis, emBasis, exBasis, alpha, beta, beta, eta, 'maxIter',2500,'rescaleRho',false);


    measValsEst = reflValsEst + flValsEst + cameraOffset;

    [pixelErr(i), pixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals,[nChannels*nFilters,nSamples]), 'absolute');

    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, reflRef, 'absolute');

    [dMatErr(i), dMatStd(i)] = fiComputeError(dMatEst, dMatRef, 'normalized');
                  

end

try
    delete(pool);
catch
end

%% Save file

if ~isempty(saveFName)
    save(saveFName,'pixelErr','pixelStd','dMatErr','dMatStd','reflErr','reflStd',...
        'exBasisGrid','emBasisGrid','alpha','beta','eta','dMatRef','reflRef',...
        'nReflBasis');
end


