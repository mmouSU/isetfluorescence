%% init
ieInit;

%%
filePath355 = which('Keratin_355nm.mat');
filePath375 = which('Keratin_375nm.mat');
filePath405 = which('Keratin_405nm.mat');
filePath435 = which('Keratin_435nm.mat');
filePath457 = which('Keratin_457nm.mat');
filePath473 = which('Keratin_473nm.mat');

%% Read spectra
wave = 365:5:705;
mag = [1, 1.4, 2.1, 3.5, 5, 5.6];
keratin355 = ieReadSpectra(filePath355, wave) / mag(1);
keratin375 = ieReadSpectra(filePath375, wave) / mag(2);
keratin405 = ieReadSpectra(filePath405, wave) / mag(3);
keratin435 = ieReadSpectra(filePath435, wave) / mag(4);
keratin457 = ieReadSpectra(filePath457, wave) / mag(5);
keratin473 = ieReadSpectra(filePath473, wave) / mag(6);

% Concatenate vector
keratin = cat(2, keratin355, keratin375, keratin405, keratin435, keratin457,...
                keratin473);
            
oldWave = [355, 375, 405, 435, 457, 473];

%% Interpolation
keratinEEM = fiEEMInterp(keratin, 'old wave', oldWave,...
                                  'new wave', wave,...
                                  'dimension', 'excitation');
%% Create fluorophore

keratinFluoSpd = keratinEEM * spd;

ieNewGraphWin; plot(wave, keratinFluoSpd);
set(gca, 'YScale', 'log'); title('Keratin Emission at 385nm')

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