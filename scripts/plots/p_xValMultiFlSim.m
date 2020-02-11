% Generate plots showing the tuning parameter cross validation results of the
% multi fluorophore algorithm. This script generates Fig. 2 from the 
% Supplemental Material.
%
% Copyright, Henryk Blasinski 2016.

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
saveDir = fullfile('~','Desktop','Figures');
% saveDir = [];
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

% Figure style parameters
lineStyle = {'rs-','gd-','bo-'};
lw = 1;
mkSz = 4;
fs = 8;
figSize = [0 0 4.5 2.25];


% Load data
inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'results','xVal',sprintf('%s_xVal_multiFl.mat',inFName));
load(fName);

nBeta = length(betaRange);
nAlpha = length(alphaRange);
nEta = length(nuRange);



% Alpha

[pixelErrPlot, minLoc] = min(reshape(permute(totalPixelErr,[2 3 1]),[nBeta*nEta, nAlpha]));
tmp = reshape(permute(totalPixelStd,[2 3 1]),[nBeta*nEta, nAlpha]);
inds = sub2ind([nBeta*nEta, nAlpha],minLoc,1:nAlpha);
pixelErrPlotStd = tmp(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reshape(permute(reflErr,[2 3 1]),[nBeta*nEta, nAlpha]));
tmp = reshape(permute(reflStd,[2 3 1]),[nBeta*nEta, nAlpha]);
inds = sub2ind([nBeta*nEta, nAlpha],minLoc,1:nAlpha);
reflErrPlotStd = tmp(inds)/sqrt(24);

[dMatNormErrPlot, minLoc] = min(reshape(permute(dMatNormErr,[2 3 1]),[nBeta*nEta, nAlpha]));
tmp = reshape(permute(dMatNormStd,[2 3 1]),[nBeta*nEta, nAlpha]);
inds = sub2ind([nBeta*nEta, nAlpha],minLoc,1:nAlpha);
dMatNormPlotStd = tmp(inds)/sqrt(24);

figure;
hold all; grid on; box on;
pl(1) = errorbar(alphaRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(alphaRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(alphaRange,dMatNormErrPlot,dMatNormPlotStd,lineStyle{3});
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
ylim([0 0.07]);
xlim([0.9*min(alphaRange) 1.1*max(alphaRange)]);
xlabel('$\alpha$','interpreter','latex');
ylabel('$\mathrm{min}_{\beta,\eta}$(RMSE)','interpreter','latex');
set(gca,'fontSize',fs);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Donaldson matrix'},'location','northwest','interpreter','latex');
if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'multiFl_xVal_alpha.eps'));
end

% Beta

[pixelErrPlot, minLoc] = min(reshape(permute(totalPixelErr,[1 3 2]),[nAlpha*nEta, nBeta]));
tmp = reshape(permute(totalPixelStd,[1 3 2]),[nAlpha*nEta, nBeta]);
inds = sub2ind([nAlpha*nEta, nBeta],minLoc,1:nBeta);
pixelErrPlotStd = tmp(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reshape(permute(reflErr,[1 3 2]),[nAlpha*nEta, nBeta]));
tmp = reshape(permute(reflStd,[1 3 2]),[nAlpha*nEta, nBeta]);
inds = sub2ind([nAlpha*nEta, nBeta],minLoc,1:nBeta);
reflErrPlotStd = tmp(inds)/sqrt(24);

[dMatNormErrPlot, minLoc] = min(reshape(permute(dMatNormErr,[1 3 2]),[nAlpha*nEta, nBeta]));
tmp = reshape(permute(dMatNormStd,[1 3 2]),[nAlpha*nEta, nBeta]);
inds = sub2ind([nAlpha*nEta, nBeta],minLoc,1:nBeta);
dMatNormPlotStd = tmp(inds)/sqrt(24);



figure;
hold all; grid on; box on;
pl(1) = errorbar(betaRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(betaRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(betaRange,dMatNormErrPlot,dMatNormPlotStd,lineStyle{3});
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
xlim([0.9*min(betaRange) 1.1*max(betaRange)]);
ylim([0 0.07]);
xlabel('$\beta$','interpreter','latex');
ylabel('$\mathrm{min}_{\alpha,\eta}$(RMSE)','interpreter','latex');
set(gca,'fontSize',fs);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Donaldson matrix'},'location','northwest','interpreter','latex');

if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'multiFl_xVal_beta.eps'));
end

% Eta

[pixelErrPlot, minLoc] = min(reshape(permute(totalPixelErr,[1 2 3]),[nAlpha*nBeta, nEta]));
tmp = reshape(permute(totalPixelStd,[1 2 3]),[nAlpha*nBeta, nEta]);
inds = sub2ind([nAlpha*nBeta, nEta],minLoc,1:nEta);
pixelErrPlotStd = tmp(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reshape(permute(reflErr,[1 2 3]),[nAlpha*nBeta, nEta]));
tmp = reshape(permute(reflStd,[1 2 3]),[nAlpha*nBeta, nEta]);
inds = sub2ind([nAlpha*nBeta, nEta],minLoc,1:nEta);
reflErrPlotStd = tmp(inds)/sqrt(24);

[dMatNormErrPlot, minLoc] = min(reshape(permute(dMatNormErr,[1 2 3]),[nAlpha*nBeta, nEta]));
tmp = reshape(permute(dMatNormStd,[1 2 3]),[nAlpha*nBeta, nEta]);
inds = sub2ind([nAlpha*nBeta, nEta],minLoc,1:nEta);
dMatNormPlotStd = tmp(inds)/sqrt(24);



figure;
hold all; grid on; box on;
pl(1) = errorbar(nuRange,pixelErrPlot,pixelErrPlotStd,lineStyle{1});
pl(2) = errorbar(nuRange,reflErrPlot,reflErrPlotStd,lineStyle{2});
pl(3) = errorbar(nuRange,dMatNormErrPlot,dMatNormPlotStd,lineStyle{3});
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
xlim([0.9*min(nuRange) 1.1*max(nuRange)]);
ylim([0 0.07]);
xlabel('$\eta$','interpreter','latex');
ylabel('$\mathrm{min}_{\alpha,\beta}$(RMSE)','interpreter','latex');
set(gca,'fontSize',fs);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Donaldson matrix'},'location','northwest','interpreter','latex');

if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'multiFl_xVal_eta.eps'));
end