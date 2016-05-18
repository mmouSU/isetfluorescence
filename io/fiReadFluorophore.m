function [ fl ] = fiReadFluorophore( fName, varargin )

% [ fl ] = fiReadFluorophore( fName,...)
%
% Read fluorophore spectral properties from a Matlab (.mat) file fName.
% Fluorescent properties must have been saved using the fiSaveFluorophore
% function to have the expected data format.
%
% Inputs:
%   fName - path to the data file
%
% Inputs (optional):
%   'wave' - a vector of waveband samples (default = 400:10:700).
%   'qe' - fluorophore quantum efficiency (default = 1).
%
% Output:
%   fl - the fiToolbox fluorophore structure
%
% Copyright, Henryk Blasinski 2016


p = inputParser;
p.addRequired('fName',@ischar);
p.addParamValue('wave',(400:10:700)',@isvector);
p.addParamValue('qe',1,@isscalar);

p.parse(fName,varargin{:});
inputs = p.Results;


data = load(inputs.fName);

fl = fluorophoreCreate('type','custom',...
                       'wave',data.wave,...
                       'name',data.name,...
                       'solvent',data.solvent,...
                       'excitation',data.excitation,...
                       'emission',data.emission);
                   
fl = fluorophoreSet(fl,'qe',inputs.qe);
if ~isempty(inputs.wave)
    fl = fluorophoreSet(fl,'wave',inputs.wave);
end



end

