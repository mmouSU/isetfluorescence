close all;
clear variables;
clc;

ieInit;

wave = 400:10:700;
nWaves = length(wave);

% Create a reflective scene
scene = sceneCreate('macbeth d65');
ieAddObject(scene);
sceneWindow();


% Create a multi-fluorophore scene
multiFlScene = fluorescentSceneCreate('height',4,'width',6,'nFluorophores',2);
multiFlScene = fluorescentSceneSet(multiFlScene,'qe',0.5);

% Create a single fluorophore scene
singleFlScene = fluorescentSceneCreate('type','onefluorophore','fluorophoreIDs',50);
singleFlScene = fluorescentSceneSet(singleFlScene,'qe',0.5);


% Combine reflectance and fluorescence
sceneMultiFl = fiSceneAddFluorescence(scene,multiFlScene);
sceneMultiFl = sceneSet(sceneMultiFl,'name','Multi-fluorophore');

sceneSingleFl = fiSceneAddFluorescence(scene,singleFlScene);
sceneSingleFl = sceneSet(sceneSingleFl,'name','Single-fluorophore');

ieAddObject(sceneMultiFl);
ieAddObject(sceneSingleFl);
sceneWindow();
