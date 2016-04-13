function [ scene ] = fiSceneAddFluorescence( scene, flScene, varargin )

p = inputParser;
p.addParamValue('replace',false,@islogical);
p.parse(varargin{:});
inputs = p.Results;

if ~strcmp(fluorescentSceneGet(flScene,'type'),'fluorescent scene')
    error('Fluorescent scene required');
end

ill = sceneGet(scene,'illuminant');

sz = sceneGet(scene,'size');
rePhotons = sceneGet(scene,'photons');

flPhotons = fluorescentSceneGet(flScene,'photons','illuminant',ill);
flPhotons = imresize(flPhotons,sz,'nearest');

if inputs.replace == false
    scene = sceneSet(scene,'photons',rePhotons + flPhotons);
else
    scene = sceneSet(scene,'photons',flPhotons);
end

end

