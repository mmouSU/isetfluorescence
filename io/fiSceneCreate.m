function scene = fiSceneCreate(fl, pattern, illuminant)
% Create a spatial pattern of fluorophore radiances in a scene structure
%
% Syntax
%    scene = fiSceneCreate(fl,pattern)
%
% Brief
%    Create a fluorescent scene with a particular spatial pattern
%
% Inputs
%    fl:       A fluorophore struct (or a cell array of such structs) 
%    pattern:  A matrix with describing the relative amount of the
%                fluorophore at  each location.
%    illuminant:  Spectral photon distribution of the illuminant
%
% Key/value pairs
%   N/A
%
% Description
%    The fluorophore and illuminant describe the spectral radiance (in
%    photons) that will be emitted.  The pattern matrix describes how much
%    of the fluorophore wil be emitted at each location.  The pattern can
%    take arbitrary positive values.
%
% Outputs
%  scene:   The scene with the fluorescent spatial pattern
%
% JEF/BW 2019
%
% See also
%

% Examples:
%{
 wave = 380:10:700;
 fName  = fullfile(fiToolboxRootPath,'data','Monici','Porphyrins.mat');
 fl  = fiReadFluorophore(fName,'wave',wave);
 % pattern = [ones(6,3) zeros(6,3)];
 slope = 2.6; imSize = 256;
 pattern = imageSlantedEdge(imSize,slope);
 ill = illuminantCreate('D65',wave);
 illuminant = illuminantGet(ill,'photons');
 scene = fiSceneCreate(fl,pattern,illuminant);
 sceneWindow(scene);
%}

%%
p = inputParser;
p.addRequired('fl',@(x)(isequal(x.type,'fluorophore')));
p.addRequired('pattern',@ismatrix);
p.addRequired('illuminant',@isvector);
p.parse(fl,pattern,illuminant);

wave = fluorophoreGet(fl,'wave');

%% Find the radiance

eem = fluorophoreGet(fl,'eem');
emission = eem*illuminant(:);
sz = size(pattern);

% Create the 3D matrix of photons, row, col, wave
photons = repmat(emission(:),1,sz(1),sz(2));
photons = permute(photons,[2 3 1]);

% Multiply by the pattern at each wavelength
for ww = 1:numel(wave)
    photons(:,:,ww) = photons(:,:,ww) .* pattern;
end

%% Put the scene together

scene = sceneCreate('empty');
scene = sceneSet(scene,'wave',wave);
scene = sceneSet(scene,'name',fluorophoreGet(fl,'name'));
scene = sceneSet(scene,'photons',photons);

end


