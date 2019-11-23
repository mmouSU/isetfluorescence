%% s_CameraFluorescence
%
% This script 
%   create a scene with reflected radiance
%   create a scene with fluorescence radiance
%   add them together 
%   simulate a camera with longpass filters and RGB sensors
% 
% SCENE
%  Simulate skin reflectance (see s_sceneReflectanceCharts.m in isetcam/scripts/scene)
%  Simulate fluorophores
%  Remember to use the same light for reflectance and fluorescence
%
% LIGHTs
%  Be sure that you use the same light for creating the radiance in the
%  reflectance scene and in the fluorescence scene to simulate a real
%  imaging scenario
%
% SENSOR
%  Things to consider when designing a camera to pick up fluorescence.
%  Shortpass filter on the light (to make sure it is very narrowband)
%  Longpass filter on the sensor (to block reflected light from the sensor)
%
% JEF, November 22 2019
%
% See also
%

%% Clean up the work space, if you like
ieInit

%% Set parameters that are shared by all scenes
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
nPixels = 64;
scene = sceneCreate('uniform equal energy',nPixels,wave);
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
fLuminance = 2;
fScene = sceneAdjustLuminance(fScene,fLuminance);
sceneWindow(fScene);

%% Combine the scenes
% Add the radiance due to reflectance and the radiance due to fluorescence
cScene = fiSceneAddFluorescence(scene, fScene );
cScene = sceneSet(cScene,'name','Combined');
sceneWindow(cScene);

%% Optics (see t_oiIntroduction)

% The simplest is method of creating an optical image is to use the
% diffraction-limited lens model.  To create a diffraction-limited
% optics with an f# of 4, you can call these functions

oi = oiCreate;
oi = oiSet(oi,'optics fnumber',4);
oi = oiSet(oi,'optics offaxis','cos4th');
oi = oiSet(oi,'optics focal length',3e-3);

% Calculate optical irradiance image
oi = oiCompute(cScene,oi);
oiWindow(oi);

%% Sensor

sensor = sensorCreate('bayer (gbrg)');  % create the sensor structure

% Set some of the key pixel properties
voltageSwing   = 1.15;  % Volts
wellCapacity   = 9000;  % Electrons
conversiongain = voltageSwing/wellCapacity;   
fillfactor     = 0.9;       % A fraction of the pixel area
pixelSize      = 2.2*1e-6;   % Meters
darkvoltage    = 1e-005;     % Volts/sec
readnoise      = 0.00096;    % Volts

% To change the fill factor, set the photodetector size and the
% pixel size to be some ratio.  To increase the fill factor 
%{
  % Compute the fill factor
  pd = sensorGet(sensor,'pixel pd size')
  pixelSize = sensorGet(sensor,'pixel size')
  fillfactor = (pd/pixelSize)^2
%}

sensorSet(sensor,'pixel pd size',[pixelSize pixelSize]*fillfactor);
sensor = sensorSet(sensor,'pixel size constant fill factor',[pixelSize pixelSize]);
sensor = sensorSet(sensor,'pixel conversion gain',conversiongain);
sensor = sensorSet(sensor,'pixel voltage swing',voltageSwing);
sensor = sensorSet(sensor,'pixel dark voltage',darkvoltage);
sensor = sensorSet(sensor,'pixel read noise volts',readnoise);

%%  Now we set some general sensor properties
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

%% Load the calibration data and attach them to the sensor structure
% Change this to be the Sony sensor
wave = sensorGet(sensor,'wave');
fullFileName = fullfile(isetRootPath,'data','sensor','colorfilters','nikon','NikonD100.mat');
[data,filterNames] = ieReadColorFilter(wave,fullFileName); 
sensor = sensorSet(sensor,'filter spectra',data);
sensor = sensorSet(sensor,'filter names',filterNames);
sensor = sensorSet(sensor,'Name','Camera-Simulation');

%% We are now ready to compute the sensor image
sensor = sensorCompute(sensor,oi);
sensorWindow(sensor);

%% Image process
ip = ipCreate;

% To see it clearly, you can do this
ip = ipCompute(ip,sensor);
ipWindow(ip);

%% To eliminate the processing you can do this
ip = ipSet(ip,'conversion method sensor','none');
ip = ipSet(ip,'internal colorspace', 'sensor');
ip = ipSet(ip,'correction method illuminant','none');
ip = ipCompute(ip,sensor);
ipWindow(ip);


%% END
