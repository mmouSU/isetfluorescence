function fl = fluorophoreSet(fl,param,val,varargin)

% Fluorophore structure setter 
%
% Copyright Henryk Blasinski, 2014

%%
if ~exist('fl','var') || isempty(fl), error('Fluorophore structure required'); end
if ~exist('param','var') || isempty(param), error('param required'); end
if ~exist('val','var') , error('val is required'); end

%%
param = ieParamFormat(param);

switch param
    case 'name'
        fl.name = val;
        
    case 'type'
        if ~strcmpi(val,'fluoropohore'), error('Type must be ''fluorophore'''); end
        fl.type = val;
        
    case 'qe'
        if (val > 1), warning('Qe greater than one, truncating to 1\n'); end
        val = min(max(val,0),1);
        
        fl.qe = val;
        
    case {'emission photons','Emission photons','emissionphotons'}
        if length(fluorophoreGet(fl,'wave')) ~= length(val), error('Wavelength sampling mismatch\n'); end
        if sum(val<0) > 0, warning('Emission less than zero, truncating\n'); end
        val = max(val,0);
        
        deltaL = fluorophoreGet(fl,'deltaWave');
        qe = 1/(sum(val)*deltaL);
        
        if qe ~= 1, warning('Emission not normalized, adjusting qe\n'); end
        fl = fluorophoreSet(fl,'qe',qe);
        
        val = val*qe;
        fl.emission = val(:);
        
    %{    
    case {'emision energy','Emission energy'}
        
        wave = fluorescenceGet(fl,'wave');
        val = Energy2Quanta(wave,val);
        
        fl = fluorescenceSet(fl,'emission photons',val);
    %}    
        
    case {'excitationphotons','excitation photons','Excitation photons'}
        
        if length(fluorophoreGet(fl,'wave')) ~= length(val), error('Wavelength sampling mismatch\n'); end
        if sum(val<0) > 0, warning('Excitation less than zero, truncating\n'); end
        val = max(val,0);
        
        if max(val) ~= 1, warning('Peak excitation different from 1, rescaling\n'); end
        val = val/max(val);
        
        fl.excitation = val(:);
        
        
      %{  
    case {'excitation energy','Excitation energy'}
        %}

    
        
    case {'wave','wavelength'}
        
        % Need to interpolate data sets and reset when wave is adjusted.
        oldW = fluorescenceGet(fl,'wave');
        newW = val(:);
        fl.spectrum.wave = newW;

        newExcitation = interp1(oldW,fluorophoreGet(fl,'excitation photons'),newW,'linear',0);
        fl = fluorophoreSet(fl,'excitation photons',newExcitation);
        
        newEmission = interp1(oldW,fluorophoreGet(fl,'emission photons'),newW,'linear',0);
        fl = fluorophoreSet(fl,'emission photons',newEmission);
        
    
    case 'solvent'
        fl.solvent = val;
        
    otherwise
        error('Unknown fluorophore parameter %s\n',param)
end

end
