function [ flSet, indx ] = fiReadFluorophoreSet( dirName, varargin )
% Read all fluorophores from a database that satisfy certain properties.
%
% Syntax:
%  [ flSet, indx ] = fiReadFluorophoreSet( dirName, ... )
% 
% Brief description:
%   Return an array of the fluorophores and their indices in the database.
%   A database index is simply the position in the sorted directory of the
%   file storing fluorescence properties.
%
% Inputs:
%   dirName - a path to the directory containing .mat files storing
%      fluorescence properties of all fluorophores.
%
% Inputs (optional):
%   'stokesShiftRange' - a (2 x 1) vector describing the minimal and
%      maximal value of the Stokes of the selected fluorophore set (default
%      = [0 Inf], i.e. all fluorophores will be returned).
%   'peakEmRange' - a (2 x 1) vector describing the minimal and maximal
%      wavelength of the emission peak (default = [0 Inf], i.e. all 
%      fluorophores will be returned).
%   'peakExRange' - a (2 x 1) vector describing the minimal and maximal
%      wavelength of the excitation peak (default = [0 Inf], i.e. all 
%      fluorophores will be returned).
%   'wave' - a (w x 1) vector of wavelength samples. All fluorescence
%      properties of all fluorophores will be re-sampled to match this
%      vector (default = 400:10:700).
%   'verbose' - a boolean indicating function verbosity. If set to true
%      certain properties of the selected fluorophores will be printed to the
%      screen (default = false).
%
% Outputs:
%   flSet - a (n x 1) array of selected fluorophore structures.
%   indx - a (n x 1) array of fluorophore indices (i.e. positions) in the
%      sorted data set directory.
%
% Copyright, Henryk Blasinski 2016.

%%
p = inputParser;
p.addRequired('dirName',@ischar);
p.addParameter('stokesShiftRange',[0 Inf],@isvector);
p.addParameter('peakEmRange',[0 Inf],@isvector);
p.addParameter('peakExRange',[0 Inf],@isvector);
p.addParameter('wave',(400:10:700)',@isvector);
p.addParameter('verbose',false,@islogical);

p.parse(dirName,varargin{:});
inputs = p.Results;

%%
files = dir(fullfile(inputs.dirName,'*.mat'));
nFiles = length(files);

flSet = [];
indx = [];
for i=1:nFiles
    
    fName = fullfile(dirName,files(i).name);
    fl = fiReadFluorophore(fName,'wave',inputs.wave);
    
    % Remove the fluorphores outside the Stokes shift range
    if fluorophoreGet(fl,'Stokes shift') < min(inputs.stokesShiftRange) ||  fluorophoreGet(fl,'Stokes shift') > max(inputs.stokesShiftRange)
        continue;
    end
    
    % Remove the fluorphores outside the emission range
    if fluorophoreGet(fl,'Peak emission') < min(inputs.peakEmRange) ||  fluorophoreGet(fl,'Peak emission') > max(inputs.peakEmRange)
        continue;
    end
    
    % Remove the fluorphores outside the emission range
    if fluorophoreGet(fl,'Peak excitation') < min(inputs.peakExRange) ||  fluorophoreGet(fl,'Peak excitation') > max(inputs.peakExRange)
        continue;
    end
    
    fl = fluorophoreSet(fl,'wave',inputs.wave);
    
    if inputs.verbose
        fprintf('%s : peak excitation %i, peak emission %i, Stokes shift %i\n',fluorophoreGet(fl,'name'),...
                                                                               fluorophoreGet(fl,'Peak excitation'),...
                                                                               fluorophoreGet(fl,'Peak emission'),...
                                                                               fluorophoreGet(fl,'Stokes shift'));
    end
    
    flSet = [flSet; fl];
    indx = [indx; i];
    
end


end

