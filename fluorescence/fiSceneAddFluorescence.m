function scene = fiSceneAddFluorescence(scene, flScene, varargin )
%
% Syntax:
%   scene = fiSceneAddFluorescence( scene, flScene, ... )
%
% Modify an Image Systems Engineering Toolbox (ISET) scene to include the
% contributions of fluorescence, as defined in the fluorescent scene. The
% illuminant used in fluorescence emission computations is defined in the
% scene structure. Note that ISET scenes are not designed to handle
% fluorescent surfaces, therefore many of the ISET functions applied to the
% returned scene structure may give incorrect results.
%
% Inputs:
%   scene - an ISET scene structure.
%   flScene - a fiToolbox fluorescent scene structure.
%
% Optional:
%   'replace' - a boolean variable indicating if the returned ISET scene
%      should contain radiance due to fluorescent emission only (default =
%      false).
%
% Returns:
%    scene - an ISET scene structure with the 'photons' field changed to
%      represent the contributions of reflected and fluoresent components
%      (if replace = false) or fluorescent component only (replace = true).
%
% Copyright, Henryk Blasinski 2016
%
% See also
%    

% Examples:
%{
scene = sceneCreate;
fScene = fiSceneCreate;
%}

%%
p = inputParser;

p.addRequired('scene', @(x)(isequal(x.type,'scene')));
p.addRequired('flScene',@(x)(isequal(x.type,'scene')));

p.addParameter('replace',false,@islogical);

p.parse(scene,flScene,varargin{:});

inputs = p.Results;

%% Get the illuminant.
ill = sceneGet(scene,'illuminant');

sz = sceneGet(scene,'size');

% Get the reflectance scene photons
rePhotons = sceneGet(scene,'photons');

% Compute the fluorescent radiance under the illuminant.
flPhotons = fluorescentSceneGet(flScene,'photons','illuminant',ill);

% Re-scale the fluorescent radiance to match the ISET scene size.
flPhotons = imresize(flPhotons,sz,'nearest');

% We either put in the sum or we just replace the reflectance with the
% fluorescence.
if inputs.replace == false
    scene = sceneSet(scene,'photons',rePhotons + flPhotons);
else
    scene = sceneSet(scene,'photons',flPhotons);
end

end

