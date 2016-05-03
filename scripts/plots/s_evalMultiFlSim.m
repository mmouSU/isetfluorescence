%% Plot the results
% Sort the fluorophores in the increasing pixel prediction error order.
close all;
clear all;
clc;

fName = fullfile(fiToolboxRootPath,'results','evaluation','McNamara-Boswell_sim_multiFl.mat');
load(fName);

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');
lw = 2;
mkSz = 8;
fs = 8;
figSize = [0 0 4.5 2.75];

% Pixel values

[sortedTotalPixelErr, indx] = sort(totalPixelErr,'ascend');
nCompounds = length(fluorophores);


figure;
hold all; grid on; box on;
pl(1) = plot(1:nCompounds,sortedTotalPixelErr,'lineWidth',lw);
pl(2) = plot(1:nCompounds,sortedTotalPixelErr - totalPixelStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(3) = plot(1:nCompounds,sortedTotalPixelErr + totalPixelStd(indx)/sqrt(24),'lineWidth',0.5*lw);
set(pl,'color','red');
lg(1) = pl(1);

% Reflectance
[sortedReflErr, indx] = sort(reflErr,'ascend');

pl(1) = plot(1:nCompounds,sortedReflErr,'lineWidth',lw);
pl(2) = plot(1:nCompounds,sortedReflErr - reflStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(3) = plot(1:nCompounds,sortedReflErr + reflStd(indx)/sqrt(24),'lineWidth',0.5*lw);
set(pl,'color','green');
lg(2) = pl(1);


% Donaldson Matrix
[sorteddMatNormErr, indx] = sort(dMatNormErr,'ascend');

pl(1) = plot(1:nCompounds,sorteddMatNormErr,'lineWidth',lw);
pl(2) = plot(1:nCompounds,sorteddMatNormErr - dMatNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(3) = plot(1:nCompounds,sorteddMatNormErr + dMatNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);
set(pl,'color','blue');
lg(3) = pl(1);


set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
xlabel('Fluorophore order');
ylabel('RMSE');
xlim([1 nCompounds]);
set(gca,'fontSize',fs);
legend(lg,{'Pixel','Reflectance','Donaldson matrix'},'location','northwest');

print('-depsc',fullfile(saveDir,'multiFl_eval.eps'));
