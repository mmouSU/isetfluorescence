close all;
clear variables;
clc;

ieInit;

wave = 400:10:700;
nWaves = length(wave);

fl = fluorophoreCreate('',wave);

ex = rand(nWaves,1);
em = rand(nWaves,1);
qe = 0.5;

fl2 = fluorophoreCreate('custom',wave,'abc','',ex,em,qe);

scene = sceneCreate('macbeth d65');
ieAddObject(scene);
sceneWindow();

sceneFl = sceneSet(scene,'fluorophore',fl);
sceneFl = sceneSet(sceneFl,'name','With fluorophore');
ieAddObject(sceneFl);
sceneWindow();


%{
fl = fluorescenceSceneCreate

fl = fluorescenceCreate();

fl = fluorescenceSet(fl,'wave',400:1:700);

fc = fluorescenceChartCreate(5,5,5,400:1:700,illuminantCreate('d65',400:1:700),[50, 100],[],1)
%}