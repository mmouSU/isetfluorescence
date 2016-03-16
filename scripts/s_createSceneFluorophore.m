close all;
clear variables;
clc;

ieInit;

wave = 400:10:700;
nWaves = length(wave);


% Create a default fluorophore
fl = fluorophoreCreate('',wave);
fl = fluorophoreSet(fl,'qe',0.1);

light = fluorophoreGet(fl,'excitation photons')*10^17;
energy = Quanta2Energy(wave,light);


% Create a reflective scene
scene = sceneCreate('macbeth d65');
ieAddObject(scene);

sceneNb = sceneAdjustIlluminant(scene,energy,0);
sceneNb = sceneSet(sceneNb,'name','Macbeth, narrowband');
ieAddObject(sceneNb);

% Create a scene with a fluorophore
scene = sceneSet(scene,'fluorophore',fl);
scene = sceneSet(scene,'name','With fluorophore, D65');
ieAddObject(scene);

% Create a scene with a fluorophore and under a narrowband illuminant
sceneFl = sceneAdjustIlluminant(scene,energy,0);
sceneFl = sceneSet(sceneFl,'name','With fluorophore, narrowband');
ieAddObject(sceneFl);



sceneWindow();


