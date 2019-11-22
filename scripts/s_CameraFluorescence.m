
%% s_CameraFluorescence

%% Clean up the work space, if you like
ieInit

%% Create a scene that has reflectance and fluorescence
wave = 380:4:780;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Grab one fluorophore
fName  = fullfile(fiToolboxRootPath,'data','Monici','Porphyrins.mat');
fl  = fiReadFluorophore(fName,'wave',wave);

fluorophorePlot(fl,'emission');
fluorophorePlot(fl,'excitation');

% Have a look at the Donaldson matrix (excitation-emission matrix)
donaldsonM = fluorophoreGet(fl,'donaldson matrix');
fluorophorePlot(fl,'donaldson matrix');

% Read in some illuminants
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
nChannels = size(illuminant,2);

% We calculate with photons and fluorescence (not energy)
illuminantPhotons = Energy2Quanta(wave,illuminant);

% Create a simple, standard scene and use one of the illuminants
scene = sceneCreate('uniform equal energy',64,wave);

whichLight = 3;    % 3 is a very short wavelength light
scene = sceneAdjustIlluminant(scene,illuminant(:,whichLight));
sceneWindow(scene);

% Calculate the fluorescence for this illuminant
illuminant = sceneGet(scene,'illuminant photons');
fName  = fullfile(fiToolboxRootPath,'data','Monici','Porphyrins.mat');
fl  = fiReadFluorophore(fName,'wave',wave);
sz = sceneGet(scene,'size');

slope = 2.6;
pattern = imageSlantedEdge(sz-1,slope);
fScene = fiSceneCreate(fl,pattern,illuminant);
sceneWindow(fScene);

% Combine the reflectance and fluorescence scene

cScene = fiSceneAddFluorescence(scene, fScene );
sceneWindow(cScene);




%% simulate camera (see t_cameraIntroduction.m)
camera = cameraCreate;

camera = cameraCompute(camera,cScene);

% Visualizing the camera objects

% To see all of the objects and their parameters listed their windows, you
% can use 
% You can open all of them with a single call, using 
cameraWindow(camera,'all');

% If you want to just see one or another, you can invoke
cameraWindow(camera,'oi');
cameraWindow(camera,'sensor');
cameraWindow(camera,'ip');
