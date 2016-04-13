function fl = fluorophoreSet(fl,param,val,varargin)

% Fluorophore structure setter 
%
% Copyright Henryk Blasinski, 2014

%%
if ~exist('fl','var') || isempty(fl), error('Fluorophore structure required'); end
if ~exist('param','var') || isempty(param), error('param required'); end
if ~exist('val','var') , error('val is required'); end

%%

% Lower case and remove spaces
param = lower(param);
param = strrep(param,' ','');

switch param
    case 'name'
        fl.name = val;
        
    case 'type'
        if ~strcmpi(val,'fluoropohore'), error('Type must be ''fluorophore'''); end
        fl.type = val;
        
    case 'qe'
        if (val > 1), warning('Qe greater than one, truncating to 1'); end
        val = min(max(val,0),1);
        
        % We only set qe for the excitation-emission representation
        if ~isfield(fl,'donaldsonMatrix')
            fl.qe = val;
        end
        
    case {'emission photons','Emission photons','emissionphotons'}
        if length(fluorophoreGet(fl,'wave')) ~= length(val), error('Wavelength sampling mismatch'); end
        
        % If the fluorophore happened to be defined with a Donaldson
        % matrix, remove the matrix from the structure;
        if isfield(fl,'donaldsonMatrix'), 
            fl = rmfield(fl,'donaldsonMatrix');
        end
        
        if sum(val<0) > 0, warning('Emission less than zero, truncating'); end
        val = max(val,0);
        
        deltaL = fluorophoreGet(fl,'deltaWave');
        qe = 1/(sum(val)*deltaL);
        if qe ~= 1, warning('Emission not normalized'); end
        
        val = val*qe;
        fl.emission = val(:);
        
      
        
    case {'excitationphotons','excitation photons','Excitation photons'}
        
        if length(fluorophoreGet(fl,'wave')) ~= length(val), error('Wavelength sampling mismatch'); end
        
        % If the fluorophore happened to be defined with a Donaldson
        % matrix, remove the matrix from the structure;
        if isfield(fl,'donaldsonMatrix'), 
            fl = rmfield(fl,'donaldsonMatrix');
        end
        
        
        
        if sum(val<0) > 0, warning('Excitation less than zero, truncating'); end
        val = max(val,0);
        
        if max(val) ~= 1, warning('Peak excitation different from 1, rescaling'); end
        val = val/max(val);
        
        fl.excitation = val(:);
        
    case {'donaldsonmatrix'}
        if length(fluorophoreGet(fl,'wave')) ~= size(val,1) || length(fluorophoreGet(fl,'wave')) ~= size(val,2)
            error('Wavelength sampling mismatch'); 
        end
        
        % Remove all fields that are relevant to one
        if isfield(fl,'excitation'), 
            fl = rmfield(fl,'excitation');
            fl = rmfield(fl,'emission');
            fl = rmfield(fl,'qe');
        end
        
        fl.donaldsonMatrix = val;
        
        
    case {'wave','wavelength'}
        
        % Need to interpolate data sets and reset when wave is adjusted.
        oldW = fluorophoreGet(fl,'wave');
        newW = val(:);
        fl.spectrum.wave = newW;

        % Interpolate excitation and emission spectra or the Donaldson
        % matrix
        if isfield(fl,'donaldsonMatrix'),
            [oldWem, oldWex] = meshgrid(oldW,oldW);
            [newWem, newWex] = meshgrid(newW,newW);
            
            newDM = interp2(oldWem,oldWex,fluorophoreGet(fl,'Donaldson matrix'),newWem,newWex,'linear',0);
            fl = fluorophoreSet(fl,'Donaldson matrix',newDM);
            
        else
            newExcitation = interp1(oldW,fluorophoreGet(fl,'excitation photons'),newW,'linear',0);
            fl = fluorophoreSet(fl,'excitation photons',newExcitation);
        
            newEmission = interp1(oldW,fluorophoreGet(fl,'emission photons'),newW,'linear',0);
            fl = fluorophoreSet(fl,'emission photons',newEmission);
        end
    
    case 'solvent'
        fl.solvent = val;
        
    otherwise
        error('Unknown fluorophore parameter %s\n',param)
end

end
