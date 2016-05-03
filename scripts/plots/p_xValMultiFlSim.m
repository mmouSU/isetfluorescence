close all;
clear all;
clc;

%% Save results
inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'results','xVal',[inFName '_xVal_multiFl.mat']);
load(fName);

nBeta = length(betaRange);
nAlpha = length(alphaRange);
nNu = length(nuRange);

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');
lineStyle = {'rs-','gd-','bo-'};
lw = 2;
mkSz = 8;
fs = 8;
figSize = [0 0 4.5 2.75];

% Alpha

[pixelErrPlot, minLoc] = min(reshape(permute(totalPixelErr,[2 3 1]),[nBeta*nNu, nAlpha]));
tmp = reshape(permute(totalPixelStd,[2 3 1]),[nBeta*nNu, nAlpha]);
inds = sub2ind([nBeta*nNu, nAlpha],minLoc,1:nAlpha);
pixelErrPlotStd = tmp(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reshape(permute(reflErr,[2 3 1]),[nBeta*nNu, nAlpha]));
tmp = reshape(permute(reflStd,[2 3 1]),[nBeta*nNu, nAlpha]);
inds = sub2ind([nBeta*nNu, nAlpha],minLoc,1:nAlpha);
reflErrPlotStd = tmp(inds)/sqrt(24);

[dMatNormErrPlot, minLoc] = min(reshape(permute(dMatNormErr,[2 3 1]),[nBeta*nNu, nAlpha]));
tmp = reshape(permute(dMatNormStd,[2 3 1]),[nBeta*nNu, nAlpha]);
inds = sub2ind([nBeta*nNu, nAlpha],minLoc,1:nAlpha);
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
xlabel('\alpha');
ylabel('min_{\beta,\eta}(RMSE)');
set(gca,'fontSize',fs);
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Donaldson matrix'},'location','northwest');
print('-depsc',fullfile(saveDir,'multiFl_xVal_alpha.eps'));


% Beta

[pixelErrPlot, minLoc] = min(reshape(permute(totalPixelErr,[1 3 2]),[nAlpha*nNu, nBeta]));
tmp = reshape(permute(totalPixelStd,[1 3 2]),[nAlpha*nNu, nBeta]);
inds = sub2ind([nAlpha*nNu, nBeta],minLoc,1:nBeta);
pixelErrPlotStd = tmp(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reshape(permute(reflErr,[1 3 2]),[nAlpha*nNu, nBeta]));
tmp = reshape(permute(reflStd,[1 3 2]),[nAlpha*nNu, nBeta]);
inds = sub2ind([nAlpha*nNu, nBeta],minLoc,1:nBeta);
reflErrPlotStd = tmp(inds)/sqrt(24);

[dMatNormErrPlot, minLoc] = min(reshape(permute(dMatNormErr,[1 3 2]),[nAlpha*nNu, nBeta]));
tmp = reshape(permute(dMatNormStd,[1 3 2]),[nAlpha*nNu, nBeta]);
inds = sub2ind([nAlpha*nNu, nBeta],minLoc,1:nBeta);
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
xlabel('\beta');
ylabel('min_{\alpha,\eta}(RMSE)');
set(gca,'fontSize',fs);
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Donaldson matrix'},'location','northwest');
print('-depsc',fullfile(saveDir,'multiFl_xVal_beta.eps'));


% Nu

[pixelErrPlot, minLoc] = min(reshape(permute(totalPixelErr,[1 2 3]),[nAlpha*nBeta, nNu]));
tmp = reshape(permute(totalPixelStd,[1 2 3]),[nAlpha*nBeta, nNu]);
inds = sub2ind([nAlpha*nBeta, nNu],minLoc,1:nNu);
pixelErrPlotStd = tmp(inds)/sqrt(24);

[reflErrPlot, minLoc] = min(reshape(permute(reflErr,[1 2 3]),[nAlpha*nBeta, nNu]));
tmp = reshape(permute(reflStd,[1 2 3]),[nAlpha*nBeta, nNu]);
inds = sub2ind([nAlpha*nBeta, nNu],minLoc,1:nNu);
reflErrPlotStd = tmp(inds)/sqrt(24);

[dMatNormErrPlot, minLoc] = min(reshape(permute(dMatNormErr,[1 2 3]),[nAlpha*nBeta, nNu]));
tmp = reshape(permute(dMatNormStd,[1 2 3]),[nAlpha*nBeta, nNu]);
inds = sub2ind([nAlpha*nBeta, nNu],minLoc,1:nNu);
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
xlabel('\eta');
ylabel('min_{\alpha,\beta}(RMSE)');
set(gca,'fontSize',fs);
set(gca,'xscale','log');
legend({'Pixel','Reflectance','Donaldson matrix'},'location','northwest');
print('-depsc',fullfile(saveDir,'multiFl_xVal_eta.eps'));
