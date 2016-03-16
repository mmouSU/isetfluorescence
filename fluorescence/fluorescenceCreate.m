function fl = fluorescenceCreate(flName,wave, varargin)

% Creates a fluorescence structure
%
% Copyright Henryk Blasinski 2014

%% Initialize parameters
if ieNotDefined('flName'), flName = 'Default fluorescence'; end

fl.name = flName;
fl.type = 'fluorescence';
fl = initDefaultSpectrum(fl,'hyperspectral');
if exist('wave','var') && ~isempty(wave)
    fl.spectrum.wave = wave; 
else
    wave = fluorescenceGet(fl,'wave');
end

%% There is no default
% The absence of a default could be a problem.

switch ieParamFormat(flName)
    
    otherwise
        emWave = 550;
        em = exp(-(fl.spectrum.wave - emWave).^2/2/(15^2));
        em = em/max(em);
        
        deltaL = wave(2) - wave(1);
        absWave = 450;
        abs = exp(-(fl.spectrum.wave - absWave).^2/2/(15^2));
        abs = abs/sum(abs)/deltaL;
        
        fl = fluorescenceSet(fl,'absorption',abs);
        fl = fluorescenceSet(fl,'emission photons',em);
        
        nWaves = fluorescenceGet(fl,'nwave');
        gainMatrix = zeros(nWaves);
        
end


return;
