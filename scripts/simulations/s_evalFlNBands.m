% Evaluate the single fluorophore algorithm using simulated data captured
% with systems containing different numbers of illuminants and spectral
% channels.
%
% This script is computationally intense and may take a long time to run.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Provide an output file name to save the results.
% saveFName = fullfile(fiToolboxRootPath,'results','evaluation',sprintf('%s_simNBands_Fl.mat',dataset));
saveFName = [];

wave = 380:4:1000;
spSize = wave(end) - wave(1);
nWaves = length(wave);
deltaL = wave(2) - wave(1);

% Algorithm tuning parameters
alpha = 0.001;
beta = 0.001;

% Number of camera filters and illuminant channels
maxChannels = 30;
maxFilters = 30;

% Target properties
flQe = 0.5;
nFluorophores = 1;
dataset = 'McNamara-Boswell';
height = 4;
width = 6;
nSamples = height*width;

% Reflectance basis
nReflBasis = 5;
nEmBasis = 12;
nExBasis = 12;
reflBasis = fiCreateBasisSet('reflectance','wave',wave','n',nReflBasis);
exBasis = fiCreateBasisSet('excitation','wave',wave','n',nExBasis);
emBasis = fiCreateBasisSet('emission','wave',wave','n',nEmBasis);

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
emRef = fluorescentSceneGet(flScene,'Emission reference');
exRef = fluorescentSceneGet(flScene,'Excitation reference');


% Placeholder variables for error computations
[channelGrid, filterGrid] = meshgrid(1:maxChannels,1:maxFilters);

reflErr = zeros(maxChannels,maxFilters);
reflStd = zeros(maxChannels,maxFilters);

dMatErr = zeros(maxChannels,maxFilters);
dMatStd = zeros(maxChannels,maxFilters);

emErr = zeros(maxChannels,maxFilters);
emStd = zeros(maxChannels,maxFilters);

exErr = zeros(maxChannels,maxFilters);
exStd = zeros(maxChannels,maxFilters);

pixelErr = zeros(maxChannels,maxFilters);
pixelStd = zeros(maxChannels,maxFilters);


try
    cluster = parcluster('local');
    cluster.NumWorkers = min(numel(channelGrid),35);
    pool = parpool(cluster,cluster.NumWorkers);
catch
end

parfor i=1:numel(channelGrid)

    nFilters = filterGrid(i);
    nChannels = channelGrid(i);
    
    cameraGain = ones(nFilters,nChannels);
    cameraOffset = zeros(nFilters,nChannels);
    
    % Define a camera
    filterWidth = ceil(spSize/nFilters/deltaL);
    camera = zeros(nWaves,nFilters);
    for c=1:nFilters
        camera(min((c-1)*filterWidth+1:c*filterWidth,nWaves),c) = 1;
    end
    
    % Define illuminant
    illuminantWidth = ceil(spSize/nChannels/deltaL);
    illuminant = zeros(nWaves,nChannels);
    for c=1:nChannels
        illuminant(min((c-1)*illuminantWidth+1:c*illuminantWidth,nWaves),c) = 1;
    end
    

    [reflValsRef, flValsRef] = fiComputeReflFlContrib(camera,illuminant,cameraGain*deltaL,reflRef,dMatRef);
    measVals = reflValsRef + flValsRef;
    
    cameraGain = repmat(cameraGain,[1 1 nSamples]);
    cameraOffset = repmat(cameraOffset,[1 1 nSamples]);
    
    nF = max(max(measVals,[],1),[],2);
    nF = repmat(nF,[nFilters,nChannels,1]);
    measVals = measVals./nF;
    cameraGain = cameraGain./nF;
    
    [ reflEst, ~, emEst, ~, exEst, ~, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndFl( measVals, camera, cameraGain*deltaL, cameraOffset, illuminant,...
        reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter', 25 );
   

    measValsEst = reflValsEst + flValsEst + cameraOffset;

    [pixelErr(i), pixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals,[nChannels*nFilters,nSamples]), 'absolute');

    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, reflRef, 'absolute');

    [emErr(i), emStd(i)] = fiComputeError(emEst,emRef,'normalized');
    [exErr(i), exStd(i)] = fiComputeError(exEst,exRef,'normalized');    
    
    dMatEst = cell(nSamples,1);
    for j=1:nSamples
        dMatEst{j} = tril(emEst(:,j)*(exEst(:,j)'),-1);
    end
    [dMatErr(i), dMatStd(i)] = fiComputeError(dMatEst, dMatRef, 'normalized');
                  

end

if ~isempty(saveFName)
    save(saveFName,'pixelErr','pixelStd','dMatErr','dMatStd','reflErr','reflStd','exErr','exStd','emErr','emStd',...
        'filterGrid','channelGrid','alpha','beta','dMatRef','reflRef',...
        'nReflBasis','nExBasis','nEmBasis');
end

try
    delete(pool);
catch
end

