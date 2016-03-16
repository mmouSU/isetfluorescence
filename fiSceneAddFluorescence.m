function [ scene ] = fiSceneAddFluorescence( scene, flScene )

if ~strcmp(fluorescentSceneGet(flScene,'type'),'fluorescent scene')
    error('Fluorescent scene required');
end

ill = sceneGet(scene,'illuminant');

sz = sceneGet(scene,'size');
rePhotons = sceneGet(scene,'photons');

flScene = fluorescentSceneSet(flScene,'scenesize',sz);
flPhotons = fluorescentSceneGet(flScene,'photons',ill);

scene = sceneSet(scene,'photons',rePhotons + flPhotons);



end

