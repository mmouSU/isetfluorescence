close all;
clear all;
clc;

fName = fullfile(fiToolboxRootPath,'results','evaluation','SNR_multiFl.mat');
load(fName);

%% Average SNR across the patches for a given condition
SNRdB = zeros(nNoiseLevels,1);
for n=1:nNoiseLevels
    
    tmp = SNR{n};
    tmp = 10*log10(tmp);
    
    SNRdB(n) = mean(tmp(:));
end




%% Pixel error

pixelErr = zeros(nSamples,nNoiseLevels,nInstances);
for n=1:nNoiseLevels
    for i=1:nInstances
        est = reflValsEst{n,i} + flValsEst{n,i};
        for s=1:nSamples
            pixelErr(s,n,i) = fiComputeError(reshape(est(:,:,s),nFilters*nChannels,1),reshape(measValsNoise{n,i}(:,:,s),nFilters*nChannels,1),'');
        end
    end
end

avgPixelErr = mean(pixelErr,3);
avgPixelErr = mean(avgPixelErr,1);

stdPixelErr = std(pixelErr,[],3);
stdPixelErr = mean(stdPixelErr,1)/sqrt(nInstances);

%% Reflectance error

reflErr = zeros(nSamples,nNoiseLevels,nInstances);
for n=1:nNoiseLevels
    for i=1:nInstances
        est = reflEst{n,i};
        for s=1:nSamples
            reflErr(s,n,i) = fiComputeError(est(:,s),reflRef(:,s),'');
        end
    end
end

avgReflErr = mean(reflErr,3);
avgReflErr = mean(avgReflErr,1);

stdReflErr = std(reflErr,[],3);
stdReflErr = mean(stdReflErr,1)/sqrt(nInstances);

%% Donaldson matrix error

dMatErr = zeros(nSamples,nNoiseLevels,nInstances);
for n=1:nNoiseLevels
    for i=1:nInstances
        est = dMatEst{n,i};
        for s=1:nSamples
            dMatErr(s,n,i) = fiComputeError(est(s),dMatRef(s),'normalized');
        end
    end
end

avgDMatErr = nanmean(dMatErr,3);
avgDMatErr = nanmean(avgDMatErr,1);

stdDMatErr = nanstd(dMatErr,[],3);
stdDMatErr = nanmean(stdDMatErr,1)/sqrt(nInstances);


%% Prepare a plot

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');

lw = 1;
fs = 8;
ms = 3;
sz = [1 1 8 5];

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


