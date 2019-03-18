function fiSaveFluorophore( fName, fl, comment )
% Save an fiFluorophore structure into a Matlab .mat file. 
%
%  Brief description
%    fiSaveFluorophore( fName, fl, ...)
%
% The fluorophore structure is under development. At the moment, it has to
% contain the following fields 
%
%    .name
%    .solvent
%    .excitation
%    .emission
%    .comment
%    .wave
%
% As we get more data, including the full Excitation-Emission matrix
% (donaldsonMatrix) from the webfluor site, we may be changing the fields.
%
% In that case excitation and emission should be empty, wave should be
% empty, and there should be
%
%     exciteW
%     emitW
%     donaldsonMatrix (length(exicteW),length(emitW))
%     
% Fluorophore model:
%
%   For much of HB's work, he assumes that there is an excitation function
%   and a single emission function.  The output is always the emission
%   function, but truncated to be longer than the excitation wavelength.
%
%   That allows him to save a fluorophore using only its emission and
%   excitation spectra.  These are normalized to unit amplitude (for
%   comparison convenience) and does not include the quantum efficiency
%   parameter, which is a scaling factor depending on many physical
%   parameters (for example concentration).
%
%   A more general model would include the full Excitation-Emission
%   Spectrum (also called the Donaldson matrix).  We have obtained full
%   Donaldson matrix data from the webfluor site.  For those datasets we
%   store the matrix in the field 'donaldsonMatrix'
%
% Inputs:
%   fName - path to where the fluorophore is to be saved.
%   fl    - the fiToolbox fluorophore structure.
%
% Inputs (optional):
%   'comment' - an optional comment string.
%
% Copyright, Henryk Blasinski 2016

%%
p = inputParser;
p.addRequired('fName',@ischar);
p.addRequired('fl',@isstruct);
p.addOptional('comment','',@ischar);

p.parse(fName,fl,comment);
inputs = p.Results;

%% We assume that the fluorophore already contains all the necessary fields.

% Why don't we just save 'fl'?
%
name       = inputs.fl.name;
solvent    = inputs.fl.solvent;
excitation = inputs.fl.excitation/max(inputs.fl.excitation);
emission   = inputs.fl.emission/max(inputs.fl.emission);
comment    = inputs.comment;
wave       = inputs.fl.spectrum.wave;

save(fName,'name','solvent','excitation','emission','comment','wave');

end

