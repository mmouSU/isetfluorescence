close all;
clear all;
clc;

dataset = 'McNamara-Boswell';
saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');

lw = 1;
fs = 8;
ms = 3;
sz = [1 1 8 5];

%% Single fluorophore

fName = fullfile(fiToolboxRootPath,'results','evaluation',sprintf('%s_simSNR_Fl.mat',dataset));
load(fName);

figure;
hold on; grid on; box on;
errorbar(SNRdB,avgPixelErr,stdPixelErr,'-rs','lineWidth',lw,'markerSize',ms);
errorbar(SNRdB,avgReflErr,stdReflErr,'-gd','lineWidth',lw,'markerSize',ms);
errorbar(SNRdB,avgExErr,stdExErr,'-bo','lineWidth',lw,'markerSize',ms);
errorbar(SNRdB,avgEmErr,stdEmErr,'-c^','lineWidth',lw,'markerSize',ms);
% set(gca,'yscale','log');
set(gca,'fontsize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
xlim([min(SNRdB) max(SNRdB)]);
xlabel('SNR, dB','fontsize',fs);
ylabel('RMSE','fontsize',fs);
legend('Pixel','Reflectance','Excitation','Emission','location','northeast');

fName = fullfile(saveDir,'FlSNR.eps');
print('-depsc',fName);

%% Multi-fluorophore

fName = fullfile(fiToolboxRootPath,'results','evaluation',sprintf('%s_simSNR_multiFl.mat',dataset));
load(fName);

figure;
hold on; grid on; box on;
errorbar(SNRdB,avgPixelErr,stdPixelErr,'-rs','lineWidth',lw,'markerSize',ms);
errorbar(SNRdB,avgReflErr,stdReflErr,'-gd','lineWidth',lw,'markerSize',ms);
errorbar(SNRdB,avgDMatErr,stdDMatErr,'-bo','lineWidth',lw,'markerSize',ms);
% set(gca,'yscale','log');
set(gca,'fontsize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
xlim([min(SNRdB) max(SNRdB)]);
xlabel('SNR, dB','fontsize',fs);
ylabel('RMSE','fontsize',fs);
legend('Pixel','Reflectance','Donaldson matrix','location','northeast');

fName = fullfile(saveDir,'multiFlSNR.eps');
print('-depsc',fName);




