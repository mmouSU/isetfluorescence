function fl = fluorescenceSet(fl,param,val,varargin)
% Fluorescence structure setter 
%
% Copyright Henryk Blasinski, 2014

%%
if ~exist('fl','var') || isempty(fl), error('fluorescence structure required'); end
if ~exist('param','var') || isempty(param), error('param required'); end
if ~exist('val','var') , error('val is required'); end

%%
param = ieParamFormat(param);

switch param
    case 'name'
        fl.name = val;
    case 'type'
        if ~strcmpi(val,'fluorescence'), error('Type must be fluorescence'); end
        fl.type = val;
    case {'emission photons','Emission photons','emissionphotons'}
        fl.emission.photons = val;
        
        mn = min(val(:));
        mx = max(val(:));
        
        fl = fluorescenceSet(fl,'emission min',mn);
        fl = fluorescenceSet(fl,'emission max',mx);
    case {'emision energy','Emission energy'}
        val = Energy2Quanta;
        
        fl = fluorescenceSet(fl,'emission photons',val);
        
    case {'emission min','emissionmin'}
        fl.emission.min = val;
    case {'emission max','emissionmax'}
        fl.emission.max = val;
        
    case {'absorption','Absorption'}
        wave = fluorescenceGet(fl,'wave');
        deltaL = wave(2) - wave(1);
        
        if (sum(val)*deltaL) ~= 1 || min(val) < 0
            % warning('Fluorescence absorption spectrum outside of range, rescaling');
            val = val - min(val);
            val = val/sum(val)/deltaL;
        end
        
        fl.absorption.data = val;
        
        mn = min(val(:));
        mx = max(val(:));
        
        fl = fluorescenceSet(fl,'absorption min',mn);
        fl = fluorescenceSet(fl,'absorption max',mx);
    case {'absorption min','absorptionmin'}
        fl.absorption.min = val;
    case {'absorption max','absorptionmax'}
        fl.absorption.max = val; 
        
    case {'wave','wavelength'}
        % il = illuminantSet(il,'wave',wave)
        % Need to interpolate data sets and reset when wave is adjusted.
        oldW = fluorescenceGet(fl,'wave');
        newW = val(:);
        fl.spectrum.wave = newW;

        newAbsorption = interp1(oldW,fluorescenceGet(fl,'absorption'),newW,'linear',0);
        fl = fluorescenceSet(fl,'absorption',newAbsorption);
        
        newEmission = interp1(oldW,fluorescenceGet(fl,'emission photons'),newW,'linear',0);
        fl = fluorescenceSet(fl,'emission photons',newEmission);
        
    case 'comment'
        fl.comment = val;
    otherwise
        error('Unknown fluorescence parameter %s\n',param)
end

end
