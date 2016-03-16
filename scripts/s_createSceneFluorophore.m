close all;
clear variables;
clc;

ieInit;

wave = 400:10:700;
nWaves = length(wave);


% Create a default fluorophore
fl = fluorophoreCreate('wave',wave);
fl = fluorophoreSet(fl,'qe',0.01);


% Create a reflective scene
scene = sceneCreate('macbeth d65');
ieAddObject(scene);

scene = fiSceneAddFluorophore(scene,fl);
ieAddObject(scene);

sceneWindow();


