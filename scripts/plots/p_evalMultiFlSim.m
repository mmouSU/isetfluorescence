% Analyze the multi fluorophore estimation results applied to the entire
% McNamara-Boswell data set. This script generates Fig. 4 in Supplemental
% Material.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
saveDir = [];
saveDir = fullfile('~','Desktop','Figures');
if ~exist(saveDir,'dir'), mkdir(saveDir); end

% Figure render parameters
lw = 1;
mkSz = 4;
fs = 8;
figSize = [0 0 4.5 2.25];

fName = fullfile(fiToolboxRootPath,'results','evaluation','McNamara-Boswell_sim_multiFl.mat');
load(fName);

%% Pixel values

[sortedTotalPixelErr, indx] = sort(totalPixelErr,'ascend');
nCompounds = length(fluorophores);


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


%% Donaldson Matrix
[sorteddMatNormErr, indx] = sort(dMatNormErr,'ascend');

pl(1) = plot(1:nCompounds,sorteddMatNormErr,'lineWidth',lw);
pl(2) = plot(1:nCompounds,sorteddMatNormErr - dMatNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(3) = plot(1:nCompounds,sorteddMatNormErr + dMatNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);
set(pl,'color','blue');
lg(3) = pl(1);


set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
set(gca,'TickLabelInterpreter','latex');
xlabel('Error rank order','interpreter','latex');
ylabel('RMSE','interpreter','latex');
xlim([1 nCompounds]);
set(gca,'fontSize',fs);
legend(lg,{'Pixel','Reflectance','Donaldson matrix'},'location','northwest','interpreter','latex');

if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'multiFl_eval.eps'));
end
