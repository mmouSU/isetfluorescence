%% Create a simple, fluorescent scene
%
%    Deprecate?
%
% Understand the basic functions and parameters
%
% JEF/BW Vistasoft, 2018
%
% See instead
%   s_sceneFluorescence, s_flSceneExample
%


%% Scene properties

dataset = 'McNamara-Boswell';
% Save results to a file saveFName
dirName = fullfile(fiToolboxRootPath,'results');
if ~exist(dirName,'dir'), mkdir(dirName); end
saveFName = fullfile(dirName,'evaluation',[dataset '_simCompare_Fl.mat']);
% saveFName = [];


% Working range for the oral eye camera
wave = 380:4:780;
deltaL = wave(2) - wave(1);
nWaves = length(wave);


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

% Here we want to read just one fluorophore instead of all 490.
% Let's figure that out
[fluorophores, ids] = fiReadFluorophoreSet(setDir,'wave',wave,...
    'peakEmRange',[wave(5) wave(end-5)],...
    'peakExRange',[wave(5) wave(end-5)]);

% This should always be 1, when we are done
nCompounds = length(ids);

%{
vcNewGraphWin;
plot(fluorophores(10).spectrum.wave,fluorophores(10).emission)
hold on;
plot(fluorophores(10).spectrum.wave,fluorophores(10).excitation)
%}

%% Create reflective scene, with no fluorescence

% We use the MCC out into the IR
scene = sceneCreate('macbethEE_IR','',wave);
scene = sceneSet(scene,'fov',5);
scene = sceneSet(scene,'distance',1);
% ieAddObject(scene); sceneWindow;

fName = fullfile(isetRootPath,'data','surfaces','macbethChart');
reflRef = ieReadSpectra(fName,wave);

flQe = 2;


%% Create fluorescent scene.  This scene has no spatial structure yet.

% We are going to put this fluorophore kind of everywhere
whichFluorophore = 50;  % We have a lot of fluorophores, so pick one
flScene = fluorescentSceneCreate('type', 'fromfluorophore',...
    'wave',wave,...
    'qe',flQe,...
    'fluorophore',fluorophores(whichFluorophore));

%{
donaldsonM = fluorescentSceneGet(flScene,'Donaldson matrix');
vcNewGraphWin;
imagesc(wave,wave,donaldsonM);
xlabel('Excitation wave (nm)'); ylabel('Emission wave (nm)');
%}

% We are making the fluorescent scene data 4x6 to match the Macbeth
% chart patches.  For such a case there are 24 Donaldson matrices, and
% 24 emission and excitation spectra.
dMatRef = fluorescentSceneGet(flScene,'Donaldson reference','sceneSize',[4 6]);
exRef   = fluorescentSceneGet(flScene,'excitation reference','sceneSize',[4 6]);
emRef   = fluorescentSceneGet(flScene,'emission reference','sceneSize',[4 6]);


%% Run ISET simulations

filterID = 1;

for ch = 1:nChannels
    % Illuminant channels
    
    fprintf('Simulating filter %i illuminant %i\n',filterID,ch);
    
    % Synthesize a fluorescent scene
    localScene = sceneAdjustIlluminant(scene,illuminant(:,ch),0);    
    localScene = fiSceneAddFluorescence(localScene,flScene);
    
    localScene = sceneSet(localScene,'name',sprintf('Illuminant channel %i',ch));
    ieAddObject(localScene); sceneWindow;

end

    
%% Not implemented yet

% Create a Point Grey Camera
[sensor, optics] = createCameraModel(filterID,'wave',wave);

%
cameraExposure = zeros(nFilters,nChannels);
cameraGain = zeros(nFilters,nChannels);
cameraOffset = zeros(nFilters,nChannels);

measVals = zeros(nFilters,nChannels,24);

%%