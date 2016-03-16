close all;
clear variables;
clc;

ieInit;

wave = 400:10:700;
nWaves = length(wave);

flScene = fluorescentSceneCreate('height',4,'width',6,'nFluorophores',2);
flScene = fluorescentSceneSet(flScene,'qe',0.5);

% Create a reflective scene
scene = sceneCreate('macbeth d65');
ieAddObject(scene);

scene = fiSceneAddFluorescence(scene,flScene);
ieAddObject(scene);

sceneWindow();

