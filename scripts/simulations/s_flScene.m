% Create a simple, fluorescent scene
%
% Understand the basic functions and parameters
%
% JEF/BW Vistasoft, 2018


%% Scene properties
dataset = 'McNamara-Boswell';
flQe = 0.1;

% Save results to a file saveFName
dirName = fullfile(fiToolboxRootPath,'results');
if ~exist(dirName,'dir'), mkdir(dirName); end
saveFName = fullfile(dirName,'evaluation',[dataset '_simCompare_Fl.mat']);
% saveFName = [];

% Working range for the oral eye camera
wave = 380:4:780;
deltaL = wave(2) - wave(1);
nWaves = length(wave);


%%
%{
% Algorithm tuning parameters
alpha = 0.1;
beta = 0.1;

alphaMultistep = 0.01;
gammaMultistep = 0.01;
%}

%% Create basis function sets
nReflBasis = 5;
nExBasis = 5;
nEmBasis = 5;

[reflBasis, reflScore] = fiCreateBasisSet('reflectance',...
    'wave',wave','n',nReflBasis);
vcNewGraphWin;
plot(wave,reflBasis,'k-');
xlabel('Wave (nm)');
grid on; 

[exBasis, exScore] = fiCreateBasisSet('excitation','wave',wave','n',nExBasis);
vcNewGraphWin;
plot(wave,exBasis,'k-');
xlabel('Wave (nm)');
grid on; 

[emBasis, emScore] = fiCreateBasisSet('emission','wave',wave','n',nEmBasis);
vcNewGraphWin;
plot(wave,emBasis,'k-');
xlabel('Wave (nm)');
grid on; 

%% Load the light spectra (in photons)
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
illuminantPhotons = Energy2Quanta(wave,illuminant);
nChannels = size(illuminant,2);

%% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);
       

%% Define which fluorophores will be used
setDir = fullfile(fiToolboxRootPath,'data',dataset);
[fluorophores, ids] = fiReadFluorophoreSet(setDir,'wave',wave,...
                           'peakEmRange',[wave(5) wave(end-5)],...
                           'peakExRange',[wave(5) wave(end-5)]);
nCompounds = length(ids);


%% Create reflective scene
scene = sceneCreate('macbethEE_IR','',wave);
scene = sceneSet(scene,'fov',5);
scene = sceneSet(scene,'distance',1);
% ieAddObject(scene); sceneWindow;

fName = fullfile(isetRootPath,'data','surfaces','macbethChart');
reflRef = ieReadSpectra(fName,wave);


%% Create a database of emission spectra
fName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell');
[flSet, fluorophoreIDs] = fiReadFluorophoreSet(fName,'wave',wave,...
            'peakEmRange',[wave(5) wave(end-5)],...
            'peakExRange',[wave(5) wave(end-5)]);

tmpScene = fluorescentSceneCreate('type','fromfluorophore',...
    'fluorophore',flSet,'wave',wave);
DB = fluorescentSceneGet(tmpScene,'emissionReference');


%% The main cross-validation loop

try
    cluster = parcluster('local');
    cluster.NumWorkers = 31;
    pool = parpool(cluster,cluster.NumWorkers);
catch 
end

flQe = 0.8;

for i=1:nCompounds

    % Create fluorescent scene
    flScene = fluorescentSceneCreate('type','fromfluorophore',...
        'wave',wave,...
        'qe',flQe,...
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
            % Illuminant channels
            
            fprintf('Simulating filter %i illuminant %i\n',f,ch);
            
            % Synthesize a fluorescent scene
            localScene = sceneAdjustIlluminant(scene,illuminant(:,ch),0);
            % ieAddObject(localScene); sceneWindow;
            
            localScene = fiSceneAddFluorescence(localScene,flScene);
            
            localScene = sceneSet(localScene,'name',sprintf('Filter %i, channel %i',f,ch));
            ieAddObject(localScene);
            % sceneWindow;
        end
    end
end 
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

    
    %% Estimation (single fluorophore)
    
    cameraGain = repmat(cameraGain,[1 1 24]);
    cameraOffset = repmat(cameraOffset,[1 1 24]);
    
    
    [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndFl( measVals, camera, cameraGain*deltaL, cameraOffset, illuminantPhotons, reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter', 50 );

    
    % Evaluation
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [singleFlTotalPixelErr(i), singleFlTotalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,24]), reshape(measVals,[nChannels*nFilters,24]), 'absolute');
    [singleFlReflPixelErr(i), singleFlReflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,24]), reshape(reflValsRef,[nChannels*nFilters,24]), 'absolute');
    [singleFlFlPixelErr(i), singleFlFlPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,24]), reshape(flValsRef,[nChannels*nFilters,24]), 'absolute');
    
    [singleFlReflErr(i), singleFlReflStd(i)] = fiComputeError(reflEst, reflRef, 'absolute');
    
    [singleFlEmErr(i), singleFlEmStd(i)] = fiComputeError(emEst, emRef, 'absolute');
    [singleFlEmNormErr(i), singleFlEmNormStd(i)] = fiComputeError(emEst, emRef, 'normalized');
    
    [singleFlExErr(i), singleFlExStd(i)] = fiComputeError(exEst, exRef, 'absolute');
    [singleFlExNormErr(i), singleFlExNormStd(i)] = fiComputeError(exEst, exRef, 'normalized');
    
    
    %% Estimation (emission only)
    [ reflEst, reflCoeffs, emEst, emCoeffs, emWghts, reflValsEst, flValsEst, hist  ] = ...
    fiRecReflAndEm( measVals, camera, cameraGain*deltaL, cameraOffset,...
    illuminantPhotons, reflBasis, emBasis, alpha, beta, 'maxIter', 50);
    
    % Evaluation
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [emTotalPixelErr(i), emTotalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,24]), reshape(measVals,[nChannels*nFilters,24]), 'absolute');
    [emReflPixelErr(i), emReflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,24]), reshape(reflValsRef,[nChannels*nFilters,24]), 'absolute');
    [emFlPixelErr(i), emFlPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,24]), reshape(flValsRef,[nChannels*nFilters,24]), 'absolute');
    
    [emReflErr(i), emReflStd(i)] = fiComputeError(reflEst, reflRef, 'absolute');
    
    [emEmErr(i), emEmStd(i)] = fiComputeError(emEst, emRef, 'absolute');
    [emEmNormErr(i), emEmNormStd(i)] = fiComputeError(emEst, emRef, 'normalized');
    
    %% Estimation (multistep Fu et al.)
    
    [ reflEst, reflCoeffs, emEst, emChromaticity, exEst, exCoeffs, reflValsEst, flValsEst, hist ] = fiRecReflAndFlMultistep( measVals,...
    camera, cameraGain*deltaL, cameraOffset, illuminantPhotons, reflBasis, DB, exBasis, alphaMultistep, gammaMultistep );


    % Evaluation
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [multistepTotalPixelErr(i), multistepTotalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,24]), reshape(measVals,[nChannels*nFilters,24]), 'absolute');
    [multistepReflPixelErr(i), multistepReflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,24]), reshape(reflValsRef,[nChannels*nFilters,24]), 'absolute');
    [multistepFlPixelErr(i), multistepFlPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,24]), reshape(flValsRef,[nChannels*nFilters,24]), 'absolute');
    
    [multistepReflErr(i), multistepReflStd(i)] = fiComputeError(reflEst, reflRef, 'absolute');
    
    [multistepEmErr(i), multistepEmStd(i)] = fiComputeError(emEst, emRef, 'absolute');
    [multistepEmNormErr(i), multistepEmNormStd(i)] = fiComputeError(emEst, emRef, 'normalized');
    
    [multistepExErr(i), multistepExStd(i)] = fiComputeError(exEst, exRef, 'absolute');
    [multistepExNormErr(i), multistepExNormStd(i)] = fiComputeError(exEst, exRef, 'normalized');
    
    
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
        'singleFlTotalPixelErr','singleFlReflPixelErr','singleFlFlPixelErr','singleFlReflErr','singleFlExErr','singleFlEmErr','singleFlExNormErr','singleFlEmNormErr',...
        'singleFlTotalPixelStd','singleFlReflPixelStd','singleFlFlPixelStd','singleFlReflStd','singleFlExStd','singleFlEmStd','singleFlExNormStd','singleFlEmNormStd',...
        'emTotalPixelErr','emReflPixelErr','emFlPixelErr','emReflErr','emEmErr','emEmNormErr',...
        'emTotalPixelStd','emReflPixelStd','emFlPixelStd','emReflStd','emEmStd','emEmNormStd',...
        'multistepTotalPixelErr','multistepReflPixelErr','multistepFlPixelErr','multistepReflErr','multistepExErr','multistepEmErr','multistepExNormErr','multistepEmNormErr',...
        'multistepTotalPixelStd','multistepReflPixelStd','multistepFlPixelStd','multistepReflStd','multistepExStd','multistepEmStd','multistepExNormStd','multistepEmNormStd');
end

%{


%% Placeholder variables for error measurements
singleFlTotalPixelErr = zeros(nCompounds,1);
singleFlReflPixelErr = zeros(nCompounds,1);
singleFlFlPixelErr = zeros(nCompounds,1);
singleFlReflErr = zeros(nCompounds,1);
singleFlExErr = zeros(nCompounds,1);
singleFlExNormErr = zeros(nCompounds,1);
singleFlEmErr = zeros(nCompounds,1);
singleFlEmNormErr = zeros(nCompounds,1);

singleFlTotalPixelStd = zeros(nCompounds,1);
singleFlReflPixelStd = zeros(nCompounds,1);
singleFlFlPixelStd = zeros(nCompounds,1);
singleFlReflStd = zeros(nCompounds,1);
singleFlExStd = zeros(nCompounds,1);
singleFlExNormStd = zeros(nCompounds,1);
singleFlEmStd = zeros(nCompounds,1);
singleFlEmNormStd = zeros(nCompounds,1);

% Placeholder variables for error measurements
emTotalPixelErr = zeros(nCompounds,1);
emReflPixelErr = zeros(nCompounds,1);
emFlPixelErr = zeros(nCompounds,1);
emReflErr = zeros(nCompounds,1);
emEmErr = zeros(nCompounds,1);
emEmNormErr = zeros(nCompounds,1);

emTotalPixelStd = zeros(nCompounds,1);
emReflPixelStd = zeros(nCompounds,1);
emFlPixelStd = zeros(nCompounds,1);
emReflStd = zeros(nCompounds,1);
emEmStd = zeros(nCompounds,1);
emEmNormStd = zeros(nCompounds,1);

% Placeholder variables for error measurements
multistepTotalPixelErr = zeros(nCompounds,1);
multistepReflPixelErr = zeros(nCompounds,1);
multistepFlPixelErr = zeros(nCompounds,1);
multistepReflErr = zeros(nCompounds,1);
multistepExErr = zeros(nCompounds,1);
multistepExNormErr = zeros(nCompounds,1);
multistepEmErr = zeros(nCompounds,1);
multistepEmNormErr = zeros(nCompounds,1);

multistepTotalPixelStd = zeros(nCompounds,1);
multistepReflPixelStd = zeros(nCompounds,1);
multistepFlPixelStd = zeros(nCompounds,1);
multistepReflStd = zeros(nCompounds,1);
multistepExStd = zeros(nCompounds,1);
multistepExNormStd = zeros(nCompounds,1);
multistepEmStd = zeros(nCompounds,1);
multistepEmNormStd = zeros(nCompounds,1);
%}