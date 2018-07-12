% Create a simple, fluorescent scene
%
% Understand the basic functions and parameters
%
% JEF/BW Vistasoft, 2018


%% Scene and fluorophore properties

wave = 380:4:780;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Grab one fluorophore
fName  = fullfile(fiToolboxRootPath,'data','LifeTechnologies','phRodoRed.mat');
fl  = fiReadFluorophore(fName,'wave',wave);

vcNewGraphWin;
semilogy(wave,fl.emission,'k-')
xlabel('Wave (nm)'); ylabel('Relative emission');

%%
donaldsonM = fluorophoreGet(fl,'donaldson matrix');
vcNewGraphWin;
imagesc(wave,wave,donaldsonM);
xlabel('Wave (nm)'); ylabel('Wave (nm)');
grid on; set(gca,'YColor',[0.8 0.8 0.8]);
set(gca,'XColor',[0.8 0.8 0.8])

%% Create a simple, standard scene
scene = sceneCreate('macbeth',[],wave);
ieAddObject(scene); sceneWindow;

%%

illuminant = sceneGet(scene,'illuminant energy');
emission = donaldsonM * illuminant(:);
vcNewGraphWin;
plot(wave,emission);

%% Make a scene that has the emission spectrum at every location
sz = sceneGet(scene,'size');

% Make an XW format of the scene energy
sceneEnergy = repmat(emission(:)',sz(1)*sz(2),1);

% Make a random amount of the fluorophore at each location
% This qe controls the spatial structure of the scene.
qe = randn(sz)*0.1 + 0.5;
qe(:,1:48) = 0;
sceneEnergy = diag(qe(:))*sceneEnergy;
sceneEnergy = XW2RGBFormat(sceneEnergy,sz(1),sz(2));

flScene = sceneCreate('macbeth',[],wave);
flScene = sceneSet(flScene,'energy',sceneEnergy);
ieAddObject(flScene); sceneWindow;

combinedScene = sceneAdd(scene,flScene);
ieAddObject(combinedScene); sceneWindow;

%%
