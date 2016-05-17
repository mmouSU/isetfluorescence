% This script demonstrates the use of functions associated with creating, and 
% operating on a single fluorophore within the fiToolbox.
%
% Copyright, Henryk Blasinski 2016.

close all;
clear all;
clc;

ieInit;

% Read a fluorophore from a database. This puts it directly into a
% fluorophore structure.
fName = fullfile(fiToolboxRootPath,'data','LifeTechnologies','AlexaFluor555');
fl = fiReadFluorophore(fName);

% Plot excitation and emission spectra
wave = fluorophoreGet(fl,'wave');
ex = fluorophoreGet(fl,'excitation');
em = fluorophoreGet(fl,'emission');

figure; 
hold on; grid on; box on;
plot(wave,[ex em],'lineWidth',2);
xlabel('Wavelength, nm');
legend({'Excitation','Emission'},'location','northeast');
title('Excitation and emission (absolute)');


% Note that you can barely see the emission spectrum. This is due to our
% convention, where the excitation spectrum is unit normalized (i.e.
% max(ex) = 1) and the quantum efficiency is entirely embedded in the
% emission spectrum. 
%
% To visualize data more easily we can request the normalized emission and
% excitation spectra. We can also limit ourselves to the visible range.

wave = 400:5:700;
fl = fiReadFluorophore(fName,'wave',wave);
exNorm = fluorophoreGet(fl,'normalized excitation');
emNorm = fluorophoreGet(fl,'normalized emission');

figure; 
hold on; grid on; box on;
plot(wave,[exNorm emNorm],'lineWidth',2);
xlabel('Wavelength, nm');
legend({'Excitation','Emission'},'location','northeast');
title('Excitation and emission (normalized)');


% We can also create a fluorophore if we know the fluorescent properties:
% excitation and emission spectra OR the Donaldson matrix. Let's try the
% latter approach here.

% Let's create some Donaldson matrix first;
DM = tril(fluorophoreGet(fl,'emission')*fluorophoreGet(fl,'excitation')',-1);

% And now we can define a new fluorophore
fl2 = fluorophoreCreate('type','fromDonaldsonMatrix','DonaldsonMatrix',DM,'wave',wave);

% Remember a fluorophore defined with a Donaldson matrix does not contain
% excitation nor emission spectra, so
ex = fluorophoreGet(fl2,'excitation');
% returns an empty matrix.

% Finally, we can have a look at the fluoresced radiance under different
% illuminants. To try out some functionalities, let's also define a new
% wavelength sampling

newWave = 400:10:700;

ill1 = illuminantCreate('D65',newWave);
ill2 = illuminantCreate('blackbody',newWave,2000);

% To get the photons we need to pass the illuminant structure (think:
% fluorescence needs to be excited with external light).
fl2 = fluorophoreSet(fl2,'wave',newWave);
photonsIll1 = fluorophoreGet(fl2,'photons',ill1);
photonsIll2 = fluorophoreGet(fl2,'photons',ill2);

figure;
hold on; grid on; box on;
plot(newWave,[photonsIll1, photonsIll2],'lineWidth',2);
xlabel('Wavelength, nm');
ylabel('Radiance, photons/nm/m^2/sr/s');
legend({'D65','2000K'},'location','northeast');
title('Radiance');






