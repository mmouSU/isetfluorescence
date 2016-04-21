close all;
clear all;
clc;

wave = 380:4:1000;

fs = 6;
sz = [1 1 5.5 3];
sz2 = [1 1 4 2];
lw = 1;

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');

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

fName = fullfile(saveDir,'qe.eps');
print('-depsc',fName);

%% Filter transmissivity

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

fName = fullfile(saveDir,'transmissivity.eps');
print('-depsc',fName);


% All filters
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

fName = fullfile(saveDir,'allFilters.eps');
print('-depsc',fName);

%% Illuminant

fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
illuminant = Energy2Quanta(wave,illuminant);

figure; box off; hold on;
plot(wave,illuminant(:,5),'LineWidth',lw,'color','cyan');
set(gca,'XTick',[400 600 800 1000]);
set(gca,'FontSize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
xlim([min(wave) max(wave)]);
xlabel('Wavelength, nm','fontsize',fs);
ylabel('q/s/sr/nm/m^2','fontsize',fs);
title('Illuminant l(\lambda)','fontsize',fs);

fName = fullfile(saveDir,'illuminant.eps');
print('-depsc',fName);


% All illuminants
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

fName = fullfile(saveDir,'allIlluminants.eps');
print('-depsc',fName);


%% Radiance

fName = fullfile(fiToolboxRootPath,'data','macbethChart');
reflectance = ieReadSpectra(fName,wave);
re = reflectance(:,3);

dirName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell');
flSet = fiReadFluorophoreSet(dirName,'wave',wave,'peakEmRange',[0 1000],'peakExRange',[400 Inf]);
fl = flSet(40);

reRad = re.*illuminant(:,5);
ill = illuminantCreate;
ill = illuminantSet(ill,'wave',wave);
ill = illuminantSet(ill,'photons',illuminant(:,5));
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

fName = fullfile(saveDir,'radiance.eps');
print('-depsc',fName);


