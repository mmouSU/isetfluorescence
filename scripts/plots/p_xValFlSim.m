% Generate plots showing the tuning parameter cross validation results of the
% single fluorophore algorithm. This script generates Fig. 3 from the 
% Supplemental Material.
%
% Copyright, Henryk Blasinski 2016.


close all;
clear all;
clc;


% Define the directory where figures will be saved. If saveDir =[], then
% figures are not saved.
saveDir = fullfile('~','Desktop','Figures');
% saveDir = [];
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

% Figure style parameters
lineStyle = {'rs-','gd-','bo-','c^-'};
lw = 1;
mkSz = 4;
fs = 8;
figSize = [0 0 4.5 2.25];


% Load data
inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'results','xVal',sprintf('%s_xVal_Fl.mat',inFName));
load(fName);

nAlpha = length(alphaRange);
nBeta = length(betaRange);




% Alpha

[pixelErrPlot, minLoc] = min(totalPixelErr,[],2);
inds = sub2ind([nAlpha, nBeta],(1:nAlpha)',minLoc);
pixelErrPlotStd = totalPixelErr(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reflErr,[],2);
inds = sub2ind([nAlpha, nBeta],(1:nAlpha)',minLoc);
reflErrPlotStd = reflErr(inds)/sqrt(24);

[exNormErrPlot, minLoc] = min(exNormErr,[],2);
inds = sub2ind([nAlpha, nBeta],(1:nAlpha)',minLoc);
exNormPlotStd = exNormErr(inds)/sqrt(24);

[emNormErrPlot, minLoc] = min(emNormErr,[],2);
inds = sub2ind([nAlpha, nBeta],(1:nAlpha)',minLoc);
emNormPlotStd = emNormErr(inds)/sqrt(24);

figure;
hold all; grid on; box on;
pl(1) = errorbar(alphaRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(alphaRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(alphaRange,exNormErrPlot,exNormPlotStd,lineStyle{3});
pl(4) = errorbar(alphaRange,emNormErrPlot,emNormPlotStd,lineStyle{4});
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
ylim([0 0.4]);
xlim([0.9*min(alphaRange) 1.1*max(alphaRange)]);
xlabel('$\alpha$','interpreter','latex');
ylabel('$\mathrm{min}_{\beta}$(RMSE)','interpreter','latex');
set(gca,'fontSize',fs);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Excitation','Emission'},'location','northwest','interpreter','latex','fontsize',fs);

if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'singleFl_xVal_alpha.eps'));
end

% Beta

[pixelErrPlot, minLoc] = min(totalPixelErr,[],1);
inds = sub2ind([nAlpha, nBeta],minLoc,1:nBeta);
pixelErrPlotStd = totalPixelErr(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reflErr,[],1);
inds = sub2ind([nAlpha, nBeta],minLoc,1:nBeta);
reflErrPlotStd = reflErr(inds)/sqrt(24);

[exNormErrPlot, minLoc] = min(exNormErr,[],1);
inds = sub2ind([nAlpha, nBeta],minLoc,1:nBeta);
exNormPlotStd = exNormErr(inds)/sqrt(24);

[emNormErrPlot, minLoc] = min(emNormErr,[],1);
inds = sub2ind([nAlpha, nBeta],minLoc,1:nBeta);
emNormPlotStd = emNormErr(inds)/sqrt(24);



figure;
hold all; grid on; box on;
pl(1) = errorbar(betaRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(betaRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(betaRange,exNormErrPlot,exNormPlotStd,lineStyle{3});
pl(4) = errorbar(betaRange,emNormErrPlot,emNormPlotStd,lineStyle{4});
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
xlim([0.9*min(betaRange) 1.1*max(betaRange)]);
ylim([0 0.4]);
xlabel('$\beta$','interpreter','latex');
ylabel('$\mathrm{min}_{\alpha}$(RMSE)','interpreter','latex');
set(gca,'fontSize',fs);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Excitation','Emission'},'location','northwest','interpreter','latex','fontsize',fs);

if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'singleFl_xVal_beta.eps'));
end


