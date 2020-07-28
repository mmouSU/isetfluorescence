function combinedScene = sceneAddFluorophore(scene,fluorophore,concentration)
% Add a fluorophore to the whole scene
%
% Synopsis
%   combinedScene = sceneAddFluorophore(scene,fluorophore,concentration)
%
% Brief description
%   For now we add it to the whole scene.  We might add it more localized
%   later.
%
% Inputs
%   scene:   Reference scene.  THe fluorescence emission is added
%            everywhere onto this scene
%  fluorophore:    See isetfluorescence
%  concetration: 
%
% Options
%  N/A
%
% Outputs
%   combinedScene:  THe original scene plus the fluorescence signal
%
% Description
%  We build the fluorescent scene as a match to the reference scene.  But
%  we change its energy to be just the fluorescence part
%
% See also
%  oeSensorChromaticities, s_MultispectralFluorescenceImaging
%

% To calculate the fluorescence we start with the scene illuminant.  This
% is generally the illuminant of the camera.  In the future we might decide
% to send in the illuminant separately from the scene illuminant.
illEnergy = sceneGet(scene,'illuminant energy');

% Now get the excitation-emission matrix for this fluorophore
eem = concentration*fluorophoreGet(fluorophore,'eem');
% ieNewGraphWin; imagesc(eem)

% Calculate the emissions for this excitation light at this fluorophore
% concentration 
emission = eem * illEnergy(:);
% ieNewGraphWin; plotRadiance(wave,emission)

% Build up the fluorescent scene from the emissions
sz = sceneGet(scene,'size');
photons = scenePhotonsFromVector(emission,sz(1),sz(2));
flScene = scene;
flScene = sceneSet(flScene,'photons',photons);
% max(emission(:))

% Combine the input scene with its fluorescent partner
combinedScene = sceneAdd(scene,flScene);

% Give it a name and return
combinedScene = sceneSet(combinedScene,'name',fluorophoreGet(fluorophore,'name'));
% sceneWindow(combinedScene);

end
