close all;
clear variables;
clc;

ieInit;

wave = 400:10:700;
nWaves = length(wave);


% Create a default fluorophore
fl = fluorophoreCreate('',wave);
fl = fluorophoreSet(fl,'qe',1);

light = fluorophoreGet(fl,'excitation photons')*10^17;
energy = Quanta2Energy(wave,light);


% Create a reflective scene
scene = sceneCreate('macbeth d65');
ieAddObject(scene);

scene = fiSceneAddFluorophore(scene,fl);
ieAddObject(scene);

sceneWindow();


