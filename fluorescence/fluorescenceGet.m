function val = fluorescenceGet(fl,param,varargin)
% Getter for the fluorescence structure
% 
% Copyright Henryk Blasinski, 2014

%% Parameter checking
if ~exist('fl','var') || isempty(fl), error('fluorescence structure required'); end
if ~exist('param','var') || isempty(param), error('param required'); end

val = [];

%% Main switch statement
param = ieParamFormat(param);
switch param
    case 'name'
        val = fl.name;
    case 'type'
        % Should always be 'illuminant'
        val = fl.type;
       
    case {'emission photons','Emission photons','emissionphotons'}
        
        if ~checkfields(fl,'emission','photons'), val = []; return; end
        val = fl.emission.photons;

    case {'absorption','Absorption'}
        
        if ~checkfields(fl,'absorption','data'), val = []; return; end
        val = fl.absorption.data;
        
    case 'wave'
        
        if isfield(fl,'spectrum'), val = fl.spectrum.wave;
        elseif ~isempty(varargin), val = sceneGet(varargin{1},'wave');
        end
        if isvector(val), val = val(:); end
        
    case 'nwave'
        % nWave = illuminantGet(il,'n wave');
        % Number of wavelength samples
        
        val = length(illuminantGet(fl,'wave'));
        
    case 'comment'
        val = fl.comment;

        
    otherwise
        error('Unknown fluorescence parameter %s\n',param)
end

end
