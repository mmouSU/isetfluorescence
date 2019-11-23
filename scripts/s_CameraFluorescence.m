
%% s_CameraFluorescence

% This script 
%   create a scene with reflected radiance
%   create a scene with fluorescence radiance
%   add them together 
%   simulate a camera with longpass filters and RGB sensors

% 
% SCENE
% Simulate skin reflectance (see s_sceneReflectanceCharts.m in isetcam/scripts/scene)
% Simulate fluorophores
% Remember to use the same light for reflectance and fluorescence

% LIGHTs
% Be sure that you use the same light for creating the radiance in the
% reflectance scene and in the fluorescence scene to simulate a real
% imaging scenario

% SENSOR
% Things to consider when designing a camera to pick up fluorescence.
% Shortpass filter on the light (to make sure it is very narrowband)
% Longpass filter on the sensor (to block reflected light from the sensor)

% 


%% Clean up the work space, if you like
ieInit
%% set some parameters that are shared by all scenes
wave = 380:4:780;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Read in some illuminants (might want to make more choices and more
% narrowband) or select from different illuminants
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave); % 14 lights
nChannels = size(illuminant,2);
% We calculate with photons and fluorescence (not energy)
illuminantPhotons = Energy2Quanta(wave,illuminant);
% Select the light
whichLight = 2;    % 2 is a very short wavelength light
sceneIlluminant = illuminant(:,whichLight)';
%% First create a scene that has radiance determined by reflectance only

% Create a simple, standard scene and use one of the illuminants
scene = sceneCreate('uniform equal energy',64,wave);

scene = sceneAdjustIlluminant(scene,sceneIlluminant);
sceneWindow(scene);

%% Second, create a scene that has radiance emitted by fluorescence
% You can create multiple fluorescence scenes to represent the effects of
% multiple fluorophores and then add them together

% Grab one fluorophore
fName  = fullfile(fiToolboxRootPath,'data','Monici','Porphyrins.mat');
fl  = fiReadFluorophore(fName,'wave',wave);

% fluorophorePlot(fl,'emission');
% fluorophorePlot(fl,'excitation');
% 
% % Have a look at the Donaldson matrix (excitation-emission matrix)
% donaldsonM = fluorophoreGet(fl,'donaldson matrix');
% fluorophorePlot(fl,'donaldson matrix');

% Calculate the fluorescence for this illuminant
illuminant = sceneGet(scene,'illuminant photons');
fName  = fullfile(fiToolboxRootPath,'data','Monici','Porphyrins.mat');
fl  = fiReadFluorophore(fName,'wave',wave);
sz = sceneGet(scene,'size');

slope = 2.6;
pattern = imageSlantedEdge(sz-1,slope);
fScene = fiSceneCreate(fl,pattern,sceneIlluminant);
sceneWindow(fScene);

%% Combine the scenes
% Add the radiance due to reflectance and the radiance due to fluorescence
cScene = fiSceneAddFluorescence(scene, fScene );
sceneWindow(cScene);

%% Optics
% see t_IntroductionOI.m
% The simplest is method of creating an optical image is to use 
% the diffraction-limited lens model.  To create a diffraction-limited
% optics with an f# of 4, you can call these functions
oi = oiCreate;
optics = oiGet(oi,'optics');           
optics = opticsSet(optics,'fnumber',4);
optics = opticsSet(optics,'offaxis','cos4th');
optics = opticsSet(optics,'focallength',3e-3);     
oi = oiSet(oi,'optics',optics);
% Calculate optical irradiance image
oi = oiCompute(scene,oi);
ieAddObject(oi); oiWindow;

%% Sensor

%% Sensor
sensor = sensorCreate('bayer (gbrg)');  % create the sensor structure
pixel =  sensorGet(sensor,'pixel'); % create the pixel structure
% Set some of the key pixel properties
voltageSwing   = 1.15;  % Volts
wellCapacity   = 9000;  % Electrons
conversiongain = voltageSwing/wellCapacity;   
fillfactor     = 0.45;       % A fraction of the pixel area
pixelSize      = 2.2*1e-6;   % Meters
darkvoltage    = 1e-005;     % Volts/sec
readnoise      = 0.00096;    % Volts
pixel = pixelSet(pixel,'size',[pixelSize pixelSize]);   
pixel = pixelSet(pixel,'conversiongain', conversiongain);        
pixel = pixelSet(pixel,'voltageswing',voltageSwing);                                             
pixel = pixelSet(pixel,'darkvoltage',darkvoltage) ;               
pixel = pixelSet(pixel,'readnoisevolts',readnoise);  

%  Now we set some general sensor properties
% exposureDuration = 0.030; % commented because we set autoexposure
dsnu =  0.0010;           % Volts (dark signal non-uniformity)
prnu = 0.2218;            % Percent (ranging between 0 and 100) photodetector response non-uniformity
analogGain   = 1;         % Used to adjust ISO speed
analogOffset = 0;         % Used to account for sensor black level
rows = 466;               % number of pixels in a row
cols = 642;               % number of pixels in a column

% Set these sensor properties
% sensor = sensorSet(sensor,'exposuretime',exposureDuration); % commented because we set autoexposure
sensorSet(sensor,'autoExposure',1);  
sensor = sensorSet(sensor,'rows',rows);
sensor = sensorSet(sensor,'cols',cols);
sensor = sensorSet(sensor,'dsnu level',dsnu);  
sensor = sensorSet(sensor,'prnu level',prnu); 
sensor = sensorSet(sensor,'analog Gain',analogGain);     
sensor = sensorSet(sensor,'analog Offset',analogOffset);   

% Stuff the pixel back into the sensor structure
sensor = sensorSet(sensor,'pixel',pixel);
sensor = pixelCenterFillPD(sensor,fillfactor);
wave = sensorGet(sensor,'wave');
% Load the calibration data and attach them to the sensor structure
% Change this to be the Sony sensor
fullFileName = fullfile(isetRootPath,'data','sensor','colorfilters','nikon','NikonD100.mat');
[data,filterNames] = ieReadColorFilter(wave,fullFileName); 
sensor = sensorSet(sensor,'filter spectra',data);
sensor = sensorSet(sensor,'filter names',filterNames);
sensor = sensorSet(sensor,'Name','Camera-Simulation');

% We are now ready to compute the sensor image
sensor = sensorCompute(sensor,oi);
ieAddObject(sensor); sensorImageWindow;
