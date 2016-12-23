% Evaluate the multi fluorophore algorithm on all fluorophores from the
% McNamara-Boswell dataset. For every fluorophore a target with 24 reflectance 
% patches is created. We also run the nuclear norm algorithm of Suo et al.
% for comparison.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Scene properties
dataset = 'McNamara-Boswell';
flQe = 0.1;

% Save results to a file saveFName (or don't if saveFName is empty)
dirName = fullfile(fiToolboxRootPath,'results');
if ~exist(dirName,'dir'), mkdir(dirName); end
saveFName = fullfile(dirName,'evaluation',[dataset '_simCompare_multiFl.mat']);
% saveFName = [];



wave = 380:4:1000;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Tuning parameters
alpha = 0.1;
beta = 0.1;
eta = 0.1;

% Tuning parameters (nuclear norm approach)
alphaNucNorm = 0.01;
sigmaNucNorm = 0.009;

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
illuminantPhotons = Energy2Quanta(wave,illuminant);
nChannels = size(illuminant,2);

% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);
       

% Define which fluorophores will be used
setDir = fullfile(fiToolboxRootPath,'data',dataset);
[fluorophores, ids] = fiReadFluorophoreSet(setDir,'wave',wave,...
                           'peakEmRange',[wave(5) wave(end-5)],...
                           'peakExRange',[wave(5) wave(end-5)]);
nCompounds = length(ids);


% Create reflective scene
scene = sceneCreate('macbethEE_IR','',wave);
scene = sceneSet(scene,'fov',5);
scene = sceneSet(scene,'distance',1);

fName = fullfile(isetRootPath,'data','surfaces','macbethChart');
reflRef = ieReadSpectra(fName,wave);

% Error placeholder variables
multiFlTotalPixelErr = zeros(nCompounds,1);
multiFlReflPixelErr = zeros(nCompounds,1);
multiFlFlPixelErr = zeros(nCompounds,1);
multiFlReflErr = zeros(nCompounds,1);
multiFldMatErr = zeros(nCompounds,1);
multiFldMatNormErr = zeros(nCompounds,1);

multiFlTotalPixelStd = zeros(nCompounds,1);
multiFlTeflPixelStd = zeros(nCompounds,1);
multiFlFlPixelStd = zeros(nCompounds,1);
multiFlReflStd = zeros(nCompounds,1);
multiFldMatStd = zeros(nCompounds,1);
multiFldMatNormStd = zeros(nCompounds,1);

nucNormTotalPixelErr = zeros(nCompounds,1);
nucNormTeflPixelErr = zeros(nCompounds,1);
nucNormFlPixelErr = zeros(nCompounds,1);
nucNormReflErr = zeros(nCompounds,1);
nucNormdMatErr = zeros(nCompounds,1);
nucNormdMatNormErr = zeros(nCompounds,1);

nucNormTotalPixelStd = zeros(nCompounds,1);
nucNormReflPixelStd = zeros(nCompounds,1);
nucNormFlPixelStd = zeros(nCompounds,1);
nucNormReflStd = zeros(nCompounds,1);
nucNormdMatStd = zeros(nCompounds,1);
nucNormdMatNormStd = zeros(nCompounds,1);

%% The main cross-validation loop


try
    cluster = parcluster('local');
    cluster.NumWorkers = 31;
    pool = parpool(cluster,cluster.NumWorkers);
catch 
end


parfor i=1:nCompounds

    % Create fluorescent scene
    flScene = fluorescentSceneCreate('type','fromfluorophore','wave',wave,'qe',flQe,...
                                     'fluorophore',fluorophores(i));
         
    dMatRef = fluorescentSceneGet(flScene,'Donaldson reference','sceneSize',[4 6]);        
    exRef = fluorescentSceneGet(flScene,'excitation reference','sceneSize',[4 6]); 
    emRef = fluorescentSceneGet(flScene,'emission reference','sceneSize',[4 6]); 
    
    
    %% Run ISET simulations
    
    cameraExposure = zeros(nFilters,nChannels);
    cameraGain = zeros(nFilters,nChannels);
    cameraOffset = zeros(nFilters,nChannels);
    
    measVals = zeros(nFilters,nChannels,24);
    
    for f=1:nFilters
        
        [sensor, optics] = createCameraModel(f,'wave',wave);
        
        for ch = 1:nChannels
            
            fprintf('Simulating filter %i channel %i\n',f,ch);
            
            % Synthesize a fluorescent scene
            localScene = sceneAdjustIlluminant(scene,illuminant(:,ch),0);
            localScene = fiSceneAddFluorescence(localScene,flScene);
            
            localScene = sceneSet(localScene,'name',sprintf('Filter %i, channel %i',f,ch));
            ieAddObject(localScene);
            
            
            % Compute the optical image
            oi = oiCreate();
            oi = oiSet(oi,'optics',optics);
            oi = oiSet(oi,'Name',sprintf('Filter %i, channel %i',f,ch));
            oi = oiCompute(localScene,oi);
            ieAddObject(oi);
            
            
            % Compute the sensor image
            FOV = [sceneGet(localScene,'fov horizontal') sceneGet(localScene,'fov vertical')];
            sensor = sensorSetSizeToFOV(sensor,FOV,localScene,oi);
            sensor = sensorSet(sensor,'Name',sprintf('Filter %i, channel %i',f,ch));
            % sensor = sensorSet(sensor,'noise flag',0);
            
            cameraExposure(f,ch) = autoExposure(oi,sensor,0.95,'luminance');
            sensor = sensorSet(sensor,'exposureTime',cameraExposure(f,ch));
            
            [cameraGain(f,ch), cameraOffset(f,ch)] = sensorGainAndOffset(localScene,oi,sensor);
            
            sensor = sensorCompute(sensor,oi);
            ieAddObject(sensor);
            
            
            % Read out pixel intensities
            sSize = sensorGet(sensor,'size');
            cornerPoints = [1 sSize(1); sSize(2) sSize(1); sSize(2) 1; 1 1];
            
            mVals = macbethSelect(sensor,0,1,cornerPoints);
            mVals = cellfun(@(x) nanmean(x)/(2^sensorGet(sensor,'nbits')),mVals);
            
            
            measVals(f,ch,:) = mVals;
            
        end
        
    end
    
    [reflValsRef, flValsRef] = fiComputeReflFlContrib(camera,illuminantPhotons,cameraGain*deltaL,reflRef,dMatRef);

    
    %% Estimation
    
    cameraGain = repmat(cameraGain,[1 1 24]);
    cameraOffset = repmat(cameraOffset,[1 1 24]);
    
    
    [ reflEst, reflCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMatEst, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndMultiFl( measVals, camera, illuminantPhotons, cameraGain*deltaL,...
        cameraOffset, reflBasis, emBasis, exBasis, alpha, beta, beta, eta, 'maxIter',250);
    
    
    % Evaluation
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [multiFlTotalPixelErr(i), multiFlTotalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,24]), reshape(measVals,[nChannels*nFilters,24]), 'absolute');
    [multiFlReflPixelErr(i), multiFlReflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,24]), reshape(reflValsRef,[nChannels*nFilters,24]), 'absolute');
    [multiFlFlPixelErr(i), multiFlFlPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,24]), reshape(flValsRef,[nChannels*nFilters,24]), 'absolute');
    
    [multiFlReflErr(i), multiFlReflStd(i)] = fiComputeError(reflEst, reflRef, 'absolute');
    
    [multiFldMatErr(i), multiFldMatStd(i)] = fiComputeError(dMatEst, dMatRef, 'absolute');
    [multiFldMatNormErr(i), multiFldMatNormStd(i)] = fiComputeError(dMatEst, dMatRef, 'normalized');
    
    
    %% Estimation (nuclear norm, Suo et al.)
    

    [ reflEst, emEst, exEst, dMatEst, reflValsEst, flValsEst, hist ] = fiRecReflAndFlNucNorm( measVals,...
    camera, cameraGain*deltaL, cameraOffset, illuminantPhotons, alphaNucNorm, sigmaNucNorm, 'maxIter', 250 );
    

    reflEst = cell2mat(cellfun(@diag,reflEst,'UniformOutput',false)');

    % Evaluation
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [nucNormTotalPixelErr(i), nucNormTotalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,24]), reshape(measVals,[nChannels*nFilters,24]), 'absolute');
    [nucNormReflPixelErr(i), nucNormReflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,24]), reshape(reflValsRef,[nChannels*nFilters,24]), 'absolute');
    [nucNormFlPixelErr(i), nucNormFlPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,24]), reshape(flValsRef,[nChannels*nFilters,24]), 'absolute');
    
    [nucNormReflErr(i), nucNormReflStd(i)] = fiComputeError(reflEst, reflRef, 'absolute');
    
    [nucNormdMatErr(i), nucNormdMatStd(i)] = fiComputeError(dMatEst, dMatRef, 'absolute');
    [nucNormdMatNormErr(i),nucNormdMatNormStd(i)] = fiComputeError(dMatEst, dMatRef, 'normalized');
    
    
    
    
    %% ISET cleanup, remove all objects.

    scenes = vcGetObjects('scene');
    vcDeleteSomeObjects('scene',1:length(scenes));
    oi = vcGetObjects('oi');
    vcDeleteSomeObjects('oi',1:length(oi));
    sensor = vcGetObjects('sensor');
    vcDeleteSomeObjects('sensor',1:length(sensor));

end

try
    delete(pool);
catch
end

%% Save results

if ~isempty(saveFName)
    save(saveFName,'fluorophores','ids','nCompounds',...
        'multiFlTotalPixelErr','multiFlReflPixelErr','multiFlFlPixelErr','multiFlReflErr','multiFldMatErr','multiFldMatNormErr',...
        'multiFlTotalPixelStd','multiFlReflPixelStd','multiFlFlPixelStd','multiFlReflStd','multiFldMatStd','multiFldMatNormStd',...
        'nucNormTotalPixelErr','nucNormReflPixelErr','nucNormFlPixelErr','nucNormReflErr','nucNormdMatErr','nucNormdMatNormErr',...
        'nucNormTotalPixelStd','nucNormReflPixelStd','nucNormFlPixelStd','nucNormReflStd','nucNormdMatStd','nucNormdMatNormStd');
end

       
