close all;
clear all;
clc;

ieInit;

% Read a fluorophore
fName = fullfile(fiToolboxRootPath,'data','LifeTechnologies','AlexaFluor555');
fl = fiReadFluorophore(fName);

% Plot excitation and emission spectra
wave = fluorophoreGet(fl,'wave');
figure; plot(wave,[fluorophoreGet(fl,'excitation') fluorophoreGet(fl,'emission')]);

% Read the fluorophore at defined wavelength sampling
fl = fiReadFluorophore(fName,'wave',400:10:700);
wave = fluorophoreGet(fl,'wave');
figure; plot(wave,[fluorophoreGet(fl,'norm excitation') fluorophoreGet(fl,'norm emission')]);

% Create a fluorophore from a Donaldson matrix
DM = fluorophoreGet(fl,'emission')*fluorophoreGet(fl,'excitation')';

fl2 = fluorophoreCreate('type','fromDonaldsonMatrix','DonaldsonMatrix',DM);

