% Analyze the single fluorophore estimation results applied to the entire
% McNamara-Boswell data set. This script generates Fig. 4 in Supplemental
% Material.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
% saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');
saveDir = [];

% Figure render parameters
lw = 2;
mkSz = 8;
fs = 8;
figSize = [0 0 4.5 2.75];

fName = fullfile(fiToolboxRootPath,'results','evaluation','McNamara-Boswell_sim_Fl.mat');
load(fName);



%% Pixel values

[sortedTotalPixelErr, indx] = sort(totalPixelErr,'ascend');

figure;
hold all; grid on; box on;
pl(1) = plot(1:nCompounds,sortedTotalPixelErr,'lineWidth',lw);
pl(2) = plot(1:nCompounds,sortedTotalPixelErr - totalPixelStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(3) = plot(1:nCompounds,sortedTotalPixelErr + totalPixelStd(indx)/sqrt(24),'lineWidth',0.5*lw);
set(pl,'color','red');
lg(1) = pl(1);


%% Reflectance
[sortedReflErr, indx] = sort(reflErr,'ascend');

pl(1) = plot(1:nCompounds,sortedReflErr,'lineWidth',lw);
pl(2) = plot(1:nCompounds,sortedReflErr - reflStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(3) = plot(1:nCompounds,sortedReflErr + reflStd(indx)/sqrt(24),'lineWidth',0.5*lw);
set(pl,'color','green');
lg(2) = pl(1);


%% Excitation
[sortedExNormErr, indx] = sort(exNormErr,'ascend');

pl(1) = plot(1:nCompounds,sortedExNormErr,'lineWidth',lw);
pl(2) = plot(1:nCompounds,sortedExNormErr - exNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(3) = plot(1:nCompounds,sortedExNormErr + exNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);
set(pl,'color','blue');
lg(3) = pl(1);

%% Emission
[sortedEmNormErr, indx] = sort(emNormErr,'ascend');

pl(1) = plot(1:nCompounds,sortedEmNormErr,'lineWidth',lw);
pl(2) = plot(1:nCompounds,sortedEmNormErr - emNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(3) = plot(1:nCompounds,sortedEmNormErr + emNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);
set(pl,'color','cyan');
lg(4) = pl(1);


set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
xlabel('Fluorophore order');
ylabel('RMSE');
xlim([1 nCompounds]);
ylim([0.0 0.5]);
set(gca,'fontSize',fs);
legend(lg,{'Pixel','Reflectance','Excitation','Emission'},'location','northwest');

if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'Fl_eval.eps'));
end
