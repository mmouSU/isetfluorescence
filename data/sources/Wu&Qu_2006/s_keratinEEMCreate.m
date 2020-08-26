%% s_keratinEMMCreate
%
% Create the EEM for keratin based on the might Wu and Gu
%
% Zheng Lyu
% Joyce Farrell
%
%


%% 
ieInit;

%% These are the emissions (energy, arbitrary units) when excited by different wavelengths

filePath355 = which('Keratin_355nm.mat');
filePath375 = which('Keratin_375nm.mat');
filePath405 = which('Keratin_405nm.mat');
filePath435 = which('Keratin_435nm.mat');
filePath457 = which('Keratin_457nm.mat');
filePath473 = which('Keratin_473nm.mat');

%% Read spectra and normalize them

wave = 365:5:705;

% These scaling factors are provided by the might Wu and Gu This make the
% relative scaling for different excitation wavelengths correct, but not
% the absolute units.
mag = [1, 1.4, 2.1, 3.5, 5, 5.6];

% The data here are in energy units.  Maybe we should convert to photons
% here?
keratin355 = ieReadSpectra(filePath355, wave) / mag(1);
keratin375 = ieReadSpectra(filePath375, wave) / mag(2);
keratin405 = ieReadSpectra(filePath405, wave) / mag(3);
keratin435 = ieReadSpectra(filePath435, wave) / mag(4);
keratin457 = ieReadSpectra(filePath457, wave) / mag(5);
keratin473 = ieReadSpectra(filePath473, wave) / mag(6);

% Make the matrix, but note that the wavelengths are not evenly spaced.
keratin = cat(2, keratin355, keratin375, keratin405, keratin435, keratin457,...
    keratin473);

% These are the wavelengths
oldWave = [355, 375, 405, 435, 457, 473];

%% Interpolate to a better wavelength sampling

keratinEEM = fiEEMInterp(keratin, 'old wave', oldWave,...
    'new wave', wave,...
    'dimension', 'excitation');

%% Create fluorophore

% If we had an spd, we could actually the expected emission for the OralEye
% camera light at 385 nm

[spd,~,comment] = ieReadSpectra('OralEye_385',wave);
keratinFluoSpd = keratinEEM * spd;

%% Plot

ieNewGraphWin; plot(wave, keratinFluoSpd);
set(gca, 'YScale', 'log'); title('Keratin Emission at 385nm')
xlabel('wavelength');
ylabel('Emission energy (normalized to 1)')

keratinFluo = fluorophoreCreate('type','fromdonaldsonmatrix',...
    'name','keratin',...
    'solvent','none', ...
    'wave', wave, ...
    'DonaldsonMatrix', keratinEEM);
%{

ieNewGraphWin;
% clim = [0 1];
% imagesc(wave, wave, eem, clim);
imagesc(wave, wave, keratinEEM);
colorbar
xticks(350:50:700)
yticks(350:50:700)
title('Keratin EEM')
    
%}
    
    
    %% Save
    savePath = fullfile(fiToolboxRootPath, 'data', 'Keratin', 'KeratinEEM.mat');
    save(savePath, 'keratinFluo');
    
    %% END