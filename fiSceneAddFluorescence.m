function [ scene ] = fiSceneAddFluorescence( scene, flScene )

if ~strcmp(fluorescentSceneGet(flScene,'type'),'fluorescent scene')
    error('Fluorescent scene required');
end

ill = sceneGet(scene,'illuminant');

sz = sceneGet(scene,'size');
rePhotons = sceneGet(scene,'photons');

flPhotons = fluorescentSceneGet(flScene,'photons','illuminant',ill);
flPhotons = imresize(flPhotons,sz,'nearest');

scene = sceneSet(scene,'photons',rePhotons + flPhotons);



end

