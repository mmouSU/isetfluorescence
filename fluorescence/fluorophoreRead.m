function [ fl ] = fluorophoreRead( fName, varargin )
% Read fluorophore spectral properties from a Matlab (.mat) file fName.
%
% [ fl ] = fluorophoreRead( fName,...)
%
% Fluorophore file should have been saved using fluorophoreSave.
%
% Inputs:
%   fName - path to the data file
%
% Inputs (optional):
%   wave - Vector of wavelengths for the fluorophore representation
%
% Output:
%   fl - the fiToolbox fluorophore structure
%
% Copyright, Henryk Blasinski 2016


% Examples
%{
   fName  = fullfile(fiToolboxRootPath,'data','LifeTechnologies','phRodoRed.mat');
   fl  = fluorophoreRead(fName);
   fluorophorePlot(fl,'donaldson matrix');
%}
%{
   fName  = fullfile(fiToolboxRootPath,'data','monici','FAD.mat');
   wave = 395:10:715;
   fl  = fluorophoreRead(fName,'wave',wave); 
%}

%%
p = inputParser;
p.addRequired('fName',@ischar);
p.addParameter('wave',[],@isvector);

p.parse(fName,varargin{:});
fname = p.Results.fName;
wave  = p.Results.wave;

%% Read in the fluorophore file
data = load(fname);

%% Create the fluorophore with the data from the file
if isfield(data, 'eem')
    if ~isempty(wave)
        data.eem = fiEEMInterp(data.eem,'old wave', data.wave,...
            'new wave', wave,...
            'dimension', 'both');
        data.wave = wave;
    end
    fl = fluorophoreCreate('type', 'fromdonaldsonmatrix',...
        'wave', data.wave,...
        'name', data.name,...
        'solvent', data.solvent,...
        'eem', data.eem);

else
    if ~isempty(wave)
        % Interpolate to these wavelengths
        data.excitation = interp1(data.wave,data.excitation,wave);
        data.emission   = interp1(data.wave,data.emission,wave);
        data.wave = wave;
    end
    
    % The emission and excitation are w.r.t. photons, not energy
    fl = fluorophoreCreate('type','custom',...
        'wave',      data.wave,...
        'name',      data.name,...
        'solvent',   data.solvent,...
        'excitation',data.excitation,...
        'emission',  data.emission);
end                
end

