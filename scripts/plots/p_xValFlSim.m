close all;
clear all;
clc;


inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'results','xVal',[inFName '_xVal_Fl.mat']);
load(fName);

nAlpha = length(alphaRange);
nBeta = length(betaRange);

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');
lineStyle = {'rs-','gd-','bo-','c^-'};
lw = 2;
mkSz = 8;
fs = 8;
figSize = [0 0 4.5 2.75];


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
xlabel('\alpha');
ylabel('min_{\beta}(RMSE)');
set(gca,'fontSize',fs);
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Excitation','Emission'},'location','northwest');
print('-depsc',fullfile(saveDir,'singleFl_xVal_alpha.eps'));


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
xlabel('\beta');
ylabel('min_{\alpha}(RMSE)');
set(gca,'fontSize',fs);
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Excitation','Emission'},'location','northwest');
print('-depsc',fullfile(saveDir,'singleFl_xVal_beta.eps'));



