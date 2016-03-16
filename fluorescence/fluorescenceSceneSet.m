function fl = fluorescenceSceneSet(fl,param,val,varargin)
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
        if ~strcmpi(val,'Fluorescence scene'), error('Type must be fluorescence'); end
        fl.type = val;
    case {'photons','Photons'}
        fl.data.photons = val;
        
        mn = min(val(:));
        mx = max(val(:));
        
        fl = fluorescenceSceneSet(fl,'photonMin',mn);
        fl = fluorescenceSceneSet(fl,'photonMax',mx);
        
    case {'photon min','photonmin'}
        fl.data.min = val;
    case {'photon max','photonmax'}
        fl.data.max = val;
    case {'excitation'}
        fl.data.excitation = val;
    case {'emission'}
        fl.data.emission = val;
    case {'patchsize'}
        fl.data.patchSize = val;
    case {'illuminant'}
        fl.illuminant = val;   
    case {'quantumefficiency'}
        fl.quantumEfficiency = val;
    case 'comment'
        fl.comment = val;
    otherwise
        error('Unknown fluorescence scene parameter %s\n',param)
end

end
