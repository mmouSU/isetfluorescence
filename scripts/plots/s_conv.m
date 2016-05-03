close all;
clear all;
clc;

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');

lw = 1;
fs = 8;
ms = 3;
sz = [1 1 8 5];


%% Multi fluorophore
fName = fullfile(fiToolboxRootPath,'results','evaluation','conv_multiFl');
data = load(fName);

pixelErr = zeros(data.maxIter,24);
dMatErr = zeros(data.maxIter,24);
reflErr = zeros(data.maxIter,24);

for i=1:24
    pixelErr(:,i) = data.hist{i}.pixelErr;
    dMatErr(:,i) = data.hist{i}.dMatErr;
    reflErr(:,i) = data.hist{i}.reflErr;
end

avgPixelErr = mean(pixelErr,2);
avgEmErr = mean(dMatErr,2);
avgReflErr = mean(reflErr,2);

stdPixelErr = std(pixelErr,0,2)/sqrt(24);
stdDMatErr = std(dMatErr,0,2)/sqrt(24);
stdReflErr = std(reflErr,0,2)/sqrt(24);


selPts = round(logspace(0,log10(data.maxIter),20));

figure; 
hold on; grid on; box on;
errorbar(selPts,avgPixelErr(selPts),stdPixelErr(selPts),'-rs','lineWidth',lw,'markerSize',ms);
errorbar(selPts,avgReflErr(selPts),stdReflErr(selPts),'-gd','lineWidth',lw,'markerSize',ms);
errorbar(selPts,avgEmErr(selPts),stdDMatErr(selPts),'-bo','lineWidth',lw,'markerSize',ms);
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'fontsize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
xlim([0 data.maxIter]);
ylim([0.009 1]);
xlabel('Iterations','fontsize',fs);
ylabel('RMSE','fontsize',fs);
legend('Pixel','Reflectance','Donaldson matrix','location','northeast');

fName = fullfile(saveDir,'multiFlConv.eps');
print('-depsc',fName);

%% Single fluorophore

fName = fullfile(fiToolboxRootPath,'results','evaluation','conv_Fl');
data = load(fName);

pixelErr = zeros(data.maxIter,24);
exErr = zeros(data.maxIter,24);
emErr = zeros(data.maxIter,24);
reflErr = zeros(data.maxIter,24);

for i=1:24
    
    if (length(data.hist{i}.pixelErr) - 1 ~= data.maxIter)
        data.hist{i}.pixelErr = [data.hist{i}.pixelErr; data.hist{i}.pixelErr(end)*ones(data.maxIter - length(data.hist{i}.pixelErr),1)];
    end
    if (length(data.hist{i}.reflErr) - 1 ~= data.maxIter)
        data.hist{i}.reflErr = [data.hist{i}.reflErr; data.hist{i}.reflErr(end)*ones(data.maxIter - length(data.hist{i}.reflErr),1)];
    end
    if (length(data.hist{i}.exErr-1) - 1 ~= data.maxIter)
        data.hist{i}.exErr = [data.hist{i}.exErr; data.hist{i}.exErr(end)*ones(data.maxIter - length(data.hist{i}.exErr),1)];
    end
    if (length(data.hist{i}.emErr-1) - 1 ~= data.maxIter)
        data.hist{i}.emErr = [data.hist{i}.emErr; data.hist{i}.emErr(end)*ones(data.maxIter - length(data.hist{i}.emErr),1)];
    end
    
    pixelErr(:,i) = data.hist{i}.pixelErr(1:data.maxIter);
    emErr(:,i) = data.hist{i}.emErr(1:data.maxIter);
    exErr(:,i) = data.hist{i}.exErr(1:data.maxIter);
    reflErr(:,i) = data.hist{i}.reflErr(1:data.maxIter);
end

avgPixelErr = mean(pixelErr,2);
avgEmErr = mean(emErr,2);
avgExErr = mean(exErr,2);
avgReflErr = mean(reflErr,2);

stdPixelErr = std(pixelErr,0,2)/sqrt(24);
stdEmErr = std(emErr,0,2)/sqrt(24);
stdExErr = std(exErr,0,2)/sqrt(24);
stdReflErr = std(reflErr,0,2)/sqrt(24);


selPts = round(logspace(0,log10(data.maxIter),20));

figure; 
hold on; grid on; box on;
errorbar(selPts,avgPixelErr(selPts),stdPixelErr(selPts),'-rs','lineWidth',lw,'markerSize',ms);
errorbar(selPts,avgReflErr(selPts),stdReflErr(selPts),'-gd','lineWidth',lw,'markerSize',ms);
errorbar(selPts,avgEmErr(selPts),stdEmErr(selPts),'-c^','lineWidth',lw,'markerSize',ms);
errorbar(selPts,avgExErr(selPts),stdExErr(selPts),'-bo','lineWidth',lw,'markerSize',ms);
set(gca,'yscale','log');
set(gca,'xscale','log');
set(gca,'fontsize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
xlim([0 data.maxIter]);
ylim([0.009 1]);
xlabel('Iterations','fontsize',fs);
ylabel('RMSE','fontsize',fs);
legend('Pixel','Reflectance','Excitation','Emission','location','northeast');

fName = fullfile(saveDir,'flConv.eps');
print('-depsc',fName);

