% Create a simple, fluorescent scene
%
%   Understand the basic functions and parameters
%
%   Requires the isetFluorescence toolbox
%
% JEF/BW Vistasoft, 2018
%
% See also
%  s_flScene, s_flSceneExample

%%
ieInit

%% Scene and fluorophore properties

wave = 380:4:780;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Grab one fluorophore
fName  = fullfile(isetRootPath,'data','fluorescence','phRodoRed.mat');
fl  = fiReadFluorophore(fName,'wave',wave);

fluorophorePlot(fl,'emission');
fluorophorePlot(fl,'excitation');

%% Have a look at the Donaldson matrix (excitation-emission matrix)
donaldsonM = fluorophoreGet(fl,'donaldson matrix');
fluorophorePlot(fl,'donaldson matrix');


%% Read in some illuminants
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
nChannels = size(illuminant,2);

% We calculate with photons and fluorescence (not energy)
illuminantPhotons = Energy2Quanta(wave,illuminant);

%% Create a simple, standard scene and use one of the illuminants
scene = sceneCreate('macbeth',[],wave);

whichLight = 3;    % 3 is a very short wavelength light
scene = sceneAdjustIlluminant(scene,illuminant(:,whichLight));
sceneWindow(scene);

%% Calculate the fluorescence for this illuminant

illuminant = sceneGet(scene,'illuminant photons');
emission = donaldsonM * illuminant(:);
ieNewGraphWin;
plot(wave,emission,'k-','linewidth',1);
grid on; xlabel('Photons'); ylabel('Wavelength (nm)');

%% Make a scene that has the emission spectrum at every location

% How big is the image?
sz = sceneGet(scene,'size');

% Make an XW format of the scene energy
sceneEnergy = repmat(emission(:)',sz(1)*sz(2),1);

% Make a random amount of the fluorophore at each location
% This qe controls the spatial structure of the scene.
fLevel = randn(sz)*0.1 + 0.5;
fLevel(:,1:48) = 0;
fLevel = RGB2XWFormat(fLevel);

% Convert multiply the emission spectrum at each point by the scalar in
% fLevel, the fluorescence emission level.
sceneEnergy = diag(fLevel(:))*sceneEnergy;
sceneEnergy = XW2RGBFormat(sceneEnergy,sz(1),sz(2));

%% Make the fluorescent scene
flScene = sceneCreate('macbeth',[],wave);
flScene = sceneSet(flScene,'energy',sceneEnergy);
ieAddObject(flScene); sceneWindow;

%% Combine the original scene with its fluorescent partner
combinedScene = sceneAdd(scene,flScene);
ieAddObject(combinedScene); sceneWindow;

%%
