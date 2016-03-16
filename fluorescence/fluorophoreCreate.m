function fl = fluorophoreCreate(name, wave, varargin)

% Creates a fluorophore structure.
%
% Copyright Henryk Blasinski 2014

%% Initialize parameters
if ieNotDefined('name'), name = 'Default'; end
if ieNotDefined('wave'), wave = 400:10:700; wave=wave(:); end

fl.name = name;
fl.type = 'fluorophore';
fl = initDefaultSpectrum(fl,'custom',wave);


%% There is no default
% The absence of a default could be a problem.

switch ieParamFormat(name)
    
    case 'custom'
    
        fl = fluorophoreCreate('',wave);
        fl = fluorophoreSet(fl,'name',varargin{1});
        fl = fluorophoreSet(fl,'solvent',varargin{2});
        fl = fluorophoreSet(fl,'excitation photons',varargin{3});
        fl = fluorophoreSet(fl,'emission photons',varargin{4});
        fl = fluorophoreSet(fl,'qe',varargin{5});
    
    otherwise
        
        % Create a default, idealized fluorophore with gaussian excitation
        % and emission spectra
        
        deltaL = wave(2) - wave(1);

        
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
