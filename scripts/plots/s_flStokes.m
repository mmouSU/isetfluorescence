close all;
clear all;
clc;

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');


fName = fullfile(fiToolboxRootPath,'results','evaluation','McNamara-Boswell_simStokes_Fl.mat');
load(fName);

lw = 1;
ms = 5;
sz = [1 1 9 4];
fs = 6;

colors = [1 0 0;
          0 1 0;
          0 0 1;
          0 1 1];

figure;
hold all; grid on; box on; 
colormap(colors);
barwitherr([totalPixelStd(:), reflStd(:),exNormStd(:), emNormStd(:)]/sqrt(24),[totalPixelErr(:), reflErr(:), exNormErr(:), emNormErr(:)]);
xlabel('Stokes shift interval, nm','fontsize',fs);
ylabel('Scaled RMSE','fontsize',fs);
set(gca,'XTick',1:4);
set(gca,'XTickLabel',{'0-25','25-50','50-75','75-100'});
set(gca,'fontsize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
ylim([0. 0.16]);
lg = legend({'Pixel','Reflectance','Excitation','Emission'},'location','northeast','fontsize',fs-2);
set(lg,'location','northeast');


fName = fullfile(saveDir,'fluorophoreStokes.eps');
print('-depsc',fName);
