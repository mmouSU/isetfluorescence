function val = fluorophoreGet(fl,param,varargin)
% Getter for the fluorophore structure
% 
% Copyright Henryk Blasinski, 2014

%% Parameter checking
if ~exist('fl','var') || isempty(fl), error('Fluorophore structure required'); end
if ~exist('param','var') || isempty(param), error('param required'); end

val = [];

%% Main switch statement
param = ieParamFormat(param);

switch param
    case 'name'
        val = fl.name;

    case 'type'
        % Should always be 'fluorophore'
        val = fl.type;
       
    case {'emission photons','Emission photons','emissionphotons'}
        
        if ~checkfields(fl,'emission'), val = []; return; end
        val = fl.emission;

    case {'excitation','excitationphotons'}
        
        if ~checkfields(fl,'excitation'), val = []; return; end
        val = fl.excitation;
        
    case 'wave'
        if isfield(fl,'spectrum'), val = fl.spectrum.wave; end
        if isvector(val), val = val(:); end
        
    case {'deltawave','deltaWave'}
        wave = fluorophoreGet(fl,'wave');
        val = wave(2) - wave(1);
        
        
    case {'photons'}
        illWave  = illuminantGet(varargin{1},'wave');
        illSpd = illuminantGet(varargin{1},'photons');
        
        fl = fluorophoreSet(fl,'wave',illWave);
        
        ex = fluorophoreGet(fl,'excitation photons');
        em = fluorophoreGet(fl,'emission photons');
        qe = fluorophoreGet(fl,'qe');
        deltaL = fluorophoreGet(fl,'deltaWave');
        
        % Apply the Stoke's constraint
        DM = qe*tril(em*ex',-1)*deltaL;
        
        val = DM*illSpd;
        
    case 'nwave'
        
        val = length(fluorophoreGet(fl,'wave'));
        
    case 'comment'
        val = fl.comment;
    
    case 'qe'
        val = fl.qe;
        
    case 'solvent'
        val = fl.solvent;
        
    otherwise
        error('Unknown fluorescence parameter %s\n',param)
end

end
