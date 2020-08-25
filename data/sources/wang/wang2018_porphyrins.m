%% Convert data to fluorophore format and save
%
% The data here will be stored on Flywheel along with the notes about
% how the data were obtained and put into the fluorophore.mat file.
%
% Porphyrin
% Mucosa
%
 % deprecated - remove this .. the data are stored in energy units.
 % convert to quanta when we create the EEM
%%
iWave = 300:5:900;

%% Porphyrins from Wang 2018

fName = 'PorphyrinsEmission_Wang2018.mat';
data = load(fName);
wave = data.PorphyrinsEmissionSpectra(:,1);
emission = data.PorphyrinsEmissionSpectra(:,2);
emission = interp1(wave,emission,iWave,'pchip',0);
emission = Energy2Quanta(iWave,emission(:));
emission = ieScale(emission,1);

ieNewGraphWin;
plot(iWave,emission,'k--');
xlabel('Wave (nm)');
ylabel('Normalized emission');

%% PorphyrinsAbsorptionSpectra (excitation)

fName = 'PorphyrinsExcitation_Wang2018.mat';
data = load(fName);
wave = data.PorphyrinsAbsorptionSpectra(:,1);
excitation = data.PorphyrinsAbsorptionSpectra(:,2);
excitation = interp1(wave,excitation,iWave,'pchip',0);
excitation = Energy2Quanta(iWave,excitation(:));
excitation = ieScale(excitation,1);

ieNewGraphWin;
plot(iWave,excitation,'k--');
xlabel('Wave (nm)');
ylabel('Normalized excitation');

%%  Porphyrin Fluorescence object

porphyrin = fluorophoreCreate('name','porphyrinWang');
porphyrin = fluorophoreSet(porphyrin,'wave',iWave);
porphyrin = fluorophoreSet(porphyrin,'emission photons',emission);
porphyrin = fluorophoreSet(porphyrin,'excitation photons',excitation);
dMatrix = fluorophoreGet(porphyrin,'donaldson matrix');

ieNewGraphWin;
imagesc(iWave,iWave,dMatrix);
xlabel('Wave (nm)');

%%

fName = fullfile(fiToolboxRootPath,'data','fluorescence','porphyrin_wang2018');
comment = sprintf('Scanned from Wang, et al. SPIE Antimicrobial blue light inactivation of Neisseria gonorrhoeae 2018, Figure 2 and 4');
fluorophoreSave(fName,porphyrin,comment);

%%  Teeth - Young person



