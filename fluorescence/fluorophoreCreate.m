function fl = fluorophoreCreate(varargin)

% Creates a fluorophore structure.
%
% Copyright Henryk Blasinski 2014

p = inputParser;

p.addParamValue('type','Default',@ischar);
p.addParamValue('wave',400:10:700,@isvector);
p.addParamValue('name','',@(x)(ischar(x) || isempty(x)));
p.addParamValue('solvent','',@(x)(ischar(x) || isempty(x)));
p.addParamValue('excitation',zeros(31,1),@isnumeric);
p.addParamValue('emission',zeros(31,1),@isnumeric);
p.addParamValue('qe',1,@isscalar);
p.addParamValue('DonaldsonMatrix',[],@isnumeric);

p.parse(varargin{:});
inputs = p.Results;


fl.name = inputs.name;
fl.type = 'fluorophore';
fl = initDefaultSpectrum(fl,'custom',inputs.wave);


%% There is no default
% The absence of a default could be a problem.

type = lower(inputs.type);
type = strrep(type,' ','');

switch type
    
    case 'fromdonaldsonmatrix'
        fl = fluorophoreCreate('wave',inputs.wave);
        fl = fluorophoreSet(fl,'name',inputs.name);
        fl = fluorophoreSet(fl,'solvent',inputs.solvent);
        fl = fluorophoreSet(fl,'Donaldson matrix',inputs.DonaldsonMatrix);
        fl = fluorophoreSet(fl,'qe',1);
    
    case 'custom'
    
        fl = fluorophoreCreate('wave',inputs.wave);
        fl = fluorophoreSet(fl,'name',inputs.name);
        fl = fluorophoreSet(fl,'solvent',inputs.solvent);
        fl = fluorophoreSet(fl,'excitation photons',inputs.excitation);
        fl = fluorophoreSet(fl,'emission photons',inputs.emission);
        fl = fluorophoreSet(fl,'qe',inputs.qe);
    
    otherwise
        
        % Create a default, idealized fluorophore with gaussian excitation
        % and emission spectra
        
        deltaL = inputs.wave(2) - inputs.wave(1);

        
        emWave = 550;
        em = exp(-(fl.spectrum.wave - emWave).^2/2/(15^2));
        em = em/sum(em)/deltaL;
        
        exWave = 450;
        ex = exp(-(fl.spectrum.wave - exWave).^2/2/(15^2));
        ex = ex/max(ex);
        
        fl = fluorophoreSet(fl,'excitation photons',ex);
        fl = fluorophoreSet(fl,'emission photons',em);
        fl = fluorophoreSet(fl,'qe',1);
        fl = fluorophoreSet(fl,'solvent','');
                
end


return;
