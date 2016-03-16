function fl = fluorescenceSceneCreate(flName, wave, varargin)
%
% Creates a fluorescence scene structure
%
% Copyright Henryk Blasinski 2014

%% Initialize parameters
if ieNotDefined('flName'), flName = 'Default'; end

fl.name = flName;
fl.type = 'Fluorescence scene';
fl = initDefaultSpectrum(fl,'hyperspectral');
if exist('wave','var') && ~isempty(wave)
    fl.spectrum.wave = wave; 
else
    wave = fluorescenceSceneGet(fl,'wave');
end


return;
