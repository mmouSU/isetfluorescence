function fiSaveFluorophore( fName, fl, varargin )
% Save an fiFluorophore structure into a Matlab .mat file. 
%
% Syntax:
%    fiSaveFluorophore( fName, fl, varargin)
%
% Brief description:
%  The fluorophore structure is under development. At the moment, it has to
%  contain the following fields 
%
%    .name       - Fluorophore name
%    .comment    - Often includes source
%    .solvent    - Solvent name
%    .eem        - Excitation emission matrix
%    .excitation - Excitation curve, separable case 
%    .emission   - Emission curve, separable case
%    .qe         - Quantum efficiency of excitation (Default is 1)
%    .spectrum.wave - Wavelength samples
%
% As we get more data, including the full Excitation-Emission matrix
% (eem) also called (donaldsonMatrix) from the webfluor site, we may
% be changing the field names.
%
% When we have the eem, the excitation and emission fields will be empty,
% and we specify the excitation and emission wavelengths
%
%     
% Fluorophore model:
%   For much of HB's work, he assumes that there is an excitation function
%   and a single emission function.  The output is always the emission
%   function, but truncated so that only wavelength longer than the
%   excitation wavelength (Stokes condition).
%
%   That allows him to save a fluorophore using only its emission and
%   excitation spectra.  These are normalized to unit amplitude (for
%   comparison convenience) and does not include the quantum efficiency
%   parameter, which is a scaling factor depending on many physical
%   parameters (for example concentration).
%
%   A more general model would include the full Excitation-Emission
%   Matrix (also called the Donaldson matrix).  We have obtained full
%   eem data from the webfluor site.  For those datasets we store the
%   matrix in the field 'eem'
%
% Inputs:
%   fName - path to where the fluorophore is to be saved.
%   fl    - the fiToolbox fluorophore structure.
%
% Inputs (optional):
%   'comment' - an optional comment string.
%
% See also
%    ISETCam


%%
p = inputParser;
p.addRequired('fName',@ischar);
p.addRequired('fl',@isstruct);
p.addOptional('comment','',@ischar);

p.parse(fName,fl,varargin{:});
inputs = p.Results;

%% We assume that the fluorophore already contains all the necessary fields.

% Why don't we just save 'fl'?
%
name       = inputs.fl.name;
solvent    = inputs.fl.solvent;
if isfield(inputs.fl,'eem')
    eem        = inputs.fl.eem;    % Excitation-emission matrix
else
    eem = [];
end

comment    = inputs.comment;

% Special case for HB when there are only excitation and emission vectors
excitation = inputs.fl.excitation/max(inputs.fl.excitation);
emission   = inputs.fl.emission/max(inputs.fl.emission);
wave       = inputs.fl.wave;

save(fName,'name','solvent','excitation','emission','eem','comment','wave');

end

