% s_lights_filters_flurophores

% The purpose of this script is to evaluate how lights x filters x
% fluorophores affect the relative amplitude of tissue fluorescence

% For any given fluorophore, the relative amplitude of emission (i.e.
% fluorescence) is determined by the match between the spectra of the light
% and the excitation spectra. In other words, each fluorophore will produce
% the greatest fluorescence if it is excited by a light that is a good
% match to the excitation spectra. 
% Hence, for each fluorophore, you would want to match the light spectra to the excitation spectra
% to produce the greatest emission
% An easy example is the visible red fluorescence generated when the tongue
% is illuminated with 400 nm.  This is the best wavelength to excite
% porphyrins present in the bacteria that sits on the tongue
 
%% plot and compare the filters
    
load('Y44.mat');
ieNewGraphWin;
plot(wavelength,data); hold on;
load('uv.mat');
plot(wavelength,data);

load('HoyaK2.mat');
ieNewGraphWin;
plot(wavelength,data);

%% plot the emission of the fluorophores
theseWaves = 300:5:800;
theseFluorophores= {'NADH.mat','elastin.mat','collagen1.mat','FAD.mat','protoporphyrin.mat' };
ieNewGraphWin;
for ii=1:numel(theseFluorophores)
    fName           = fullfile(fiToolboxRootPath,'data','webfluor',theseFluorophores{ii});
    fluor(ii)       = fiReadFluorophore(fName,'wave',theseWaves); 
    plot(theseWaves,fluor(ii).emission); hold on
end
legend(theseFluorophores,'FontSize',18)
title('Emission');

%% Plot the excitation of the fluorophores
theseWaves = 300:5:800;
theseFluorophores= {'NADH.mat','elastin.mat','collagen1.mat','FAD.mat','protoporphyrin.mat' };
ieNewGraphWin;
for ii=1:numel(theseFluorophores)
    fName           = fullfile(fiToolboxRootPath,'data','webfluor',theseFluorophores{ii});
    fluor(ii)       = fiReadFluorophore(fName,'wave',theseWaves); 
    plot(theseWaves,fluor(ii).excitation); hold on
end
legend(theseFluorophores,'FontSize',18)
title('Excitation');

%% Plot the lights
load('LED405.mat');
ieNewGraphWin;
plot(wavelength,data/max(data),'b','LineWidth',2);
hold on;
load('LED425.mat');
plot(wavelength,data/max(data),'g','LineWidth',2);
load('LED450.mat');
plot(wavelength,data/max(data),'r','LineWidth',2);
legend('405 nm', '425 nm', '450 nm');
title('LED spectra');
ax = gca;
ax.FontSize=16;


%% Multiply fluorophores with 475 nm Longpass filter to see how much fluorescence gets through


