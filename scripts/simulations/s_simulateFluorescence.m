close all;
clear all;
clc;

ieInit;

dataset = 'McNamara-Boswell';
flQe = 0.1;
height = 4;
width = 6;
nFluorophores = 1;

wave = 380:4:1000;
nWaves = length(wave);
deltaL = wave(2) - wave(1);


% Load the light spectra
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
nChannels = size(illuminant,2);

% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);

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
exRef = fluorescentSceneGet(flScene,'excitation reference');
emRef = fluorescentSceneGet(flScene,'emission reference');


           

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

sceneWindow;
oiWindow;
sensorWindow;

%% Compute linear model predictions

illuminantPhotons = Energy2Quanta(wave,illuminant);
[reflValsRef, flValsRef] = fiComputeReflFlContrib(camera,illuminantPhotons,cameraGain*deltaL,reflRef,dMatRef);

predVals = reflValsRef + flValsRef + repmat(cameraOffset,[1 1 24]);


figure;
hold all; grid on; box on;
plot(predVals(:),measVals(:),'.');
xlabel('Linear model pixel intensity');
ylabel('ISET simulation pixel intensity');

%% Save data

fName = fullfile(fiToolboxRootPath,'data','simulations',sprintf('%s_%ix%ix%i_qe_%0.2f.mat',dataset,height,width,nFluorophores,flQe));

save(fName,'cameraGain','cameraOffset','measVals','reflRef','dMatRef','exRef','emRef','reflValsRef','flValsRef',...
           'wave');



