% Testing whether the EEM energy and EEM photons are properly scaled.

%% Init
ieInit;

%%
wave = 365:5:705;
flName = 'elastin_webfluor';
thisFluo = fluorophoreRead(flName,'wave',wave');
eemEnergy = fluorophoreGet(thisFluo, 'eem energy');
eemPhoton = fluorophoreGet(thisFluo, 'eem photons');

% Create a spectrum in photon
spdPhoton = ones(size(wave));
spdEnergy = Quanta2Energy(wave, spdPhoton);

% The number of photons we get based on the eemPhoton will be:
emissionPhoton = eemPhoton * spdPhoton';
emissionEnergy = eemEnergy * spdEnergy';

%{
ieNewGraphWin;
hold all;
plot(wave, emissionPhoton / max(emissionPhoton(:)));
plot(wave, emissionEnergy / max(emissionEnergy(:)));
%}

% Convert emissionPhoton to emissionEnergy
emissionEnergyFromPhoton = Quanta2Energy(wave, emissionPhoton);

%{
ieNewGraphWin;
hold all;
plot(emissionEnergy, emissionEnergyFromPhoton);
identityLine;
%}

%% END