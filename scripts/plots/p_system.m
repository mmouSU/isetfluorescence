% This script generates a variety of plots of the illuminant spectral power
% distributions and filter transmissivities used in the experimental system
% we built. Some of these plots are used in Fig. 1 and others in Fig. 9 in
% the paper.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
% saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');
saveDir = [];


wave = 380:4:1000;

% Figure display settings
fs = 6;
sz = [1 1 5.5 3];
sz2 = [1 1 4 2];
lw = 1;


%% Camera qe
fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

figure; box off; hold on;
plot(wave,qe,'LineWidth',lw,'color','black');
set(gca,'YTick',[0 0.1 0.2 0.3 0.4 0.5]);
set(gca,'XTick',[400 600 800 1000]);
set(gca,'FontSize',fs-2);
set(gcf,'PaperUnits','centimeters');
xlim([min(wave) max(wave)]);
ylim([0 0.6]);
set(gcf,'PaperPosition',sz);
xlabel('Wavelength, nm','fontsize',fs);
title('Quantum efficiency q_e(\lambda)','fontsize',fs);

if ~isempty(saveDir)
    fName = fullfile(saveDir,'qe.eps');
    print('-depsc',fName);
end

%% Filter transmissivity

% Plot one specific filter
fName = fullfile(fiToolboxRootPath,'camera','filters');
camera = ieReadSpectra(fName,wave);

figure; box off; hold on;
plot(wave,camera(:,4),'LineWidth',lw,'color',[0.5 1 0.2]);
set(gca,'YTick',0:0.2:1);
set(gca,'XTick',[400 600 800 1000]);
set(gca,'FontSize',fs-2);
set(gcf,'PaperUnits','centimeters');
xlim([min(wave) max(wave)]);
ylim([0 1.05]);
set(gcf,'PaperPosition',sz);
xlabel('Wavelength, nm','fontsize',fs);
title('Transmissvity s(\lambda)','fontsize',fs);

if ~isempty(saveDir)
    fName = fullfile(saveDir,'transmissivity.eps');
    print('-depsc',fName);
end

% Plot all filters
figure; box on; grid on; hold on;
plot(wave,diag(qe)*camera,'LineWidth',lw);
set(gca,'XTick',[400 600 800 1000]);
set(gca,'FontSize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
xlim([min(wave) max(wave)]);
ylim([0 0.6]);
xlabel('Wavelength, nm','fontsize',fs);
ylabel('Transmissivity','fontsize',fs);

if ~isempty(saveDir)
    fName = fullfile(saveDir,'allFilters.eps');
    print('-depsc',fName);
end

%% Illuminant

% Plot one specific illuminant
illID = 5;

fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
illuminant = Energy2Quanta(wave,illuminant);

figure; box off; hold on;
plot(wave,illuminant(:,illID),'LineWidth',lw,'color','cyan');
set(gca,'XTick',[400 600 800 1000]);
set(gca,'FontSize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
xlim([min(wave) max(wave)]);
xlabel('Wavelength, nm','fontsize',fs);
ylabel('q/s/sr/nm/m^2','fontsize',fs);
title('Illuminant l(\lambda)','fontsize',fs);

if ~isempty(saveDir)
    fName = fullfile(saveDir,'illuminant.eps');
    print('-depsc',fName);
end

% Plot all illuminants
figure; box on; grid on; hold on;
plot(wave,illuminant/max(illuminant(:)),'LineWidth',lw);
set(gca,'XTick',[400 600 800 1000]);
set(gca,'FontSize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
xlim([min(wave) max(wave)]);
ylim([0 1.05]);
xlabel('Wavelength, nm','fontsize',fs);
ylabel('Norm. intensity','fontsize',fs);

if ~isempty(saveDir)
    fName = fullfile(saveDir,'allIlluminants.eps');
    print('-depsc',fName);
end

%% Radiance
% Compute reflected and fluoresced radiances under selected illuminants.

% Pick a reflectance
fName = fullfile(fiToolboxRootPath,'data','macbethChart');
reflectance = ieReadSpectra(fName,wave);
re = reflectance(:,3);

% Pick a random fluorophore
dirName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell');
flSet = fiReadFluorophoreSet(dirName,'wave',wave,'peakEmRange',[0 1000],'peakExRange',[400 Inf]);
fl = flSet(40);

% Compute the radiances
reRad = re.*illuminant(:,illID);
ill = illuminantCreate;
ill = illuminantSet(ill,'wave',wave);
ill = illuminantSet(ill,'photons',illuminant(:,illID));
flRad = fluorophoreGet(fl,'photons',ill);


figure; box off; hold on;
plot(wave,reRad,'LineWidth',lw);
plot(wave,flRad,'--','LineWidth',lw,'color','magenta');
set(gca,'XTick',[400 600 800 1000]);
set(gca,'FontSize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
xlim([min(wave) max(wave)]);
xlabel('Wavelength, nm','fontsize',fs);
ylabel('q/s/sr/nm/m^2','fontsize',fs);
title({'Radiance \rho(\lambda)'},'fontsize',fs);

lg = legend({'Reflected \rho_r(\lambda)', 'Fluoresced \rho_f(\lambda)'},'fontsize',fs-2);
set(lg,'Box','off');

if ~isempty(saveDir)
    fName = fullfile(saveDir,'radiance.eps');
    print('-depsc',fName);
end

