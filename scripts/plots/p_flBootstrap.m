% Plot single fluorophore and CIM estimation results. Estimate confidence
% intervals through bootstrapping. This script generates Fig. 14 in the
% paper.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
% saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');
saveDir = [];

% Figure display properties
fs = 10;
lw = 2;
ms = 7;
sz = [1 1 3.5 2.5];

% Figure legend
leg = {'Reference','Single fl.','CIM'};

% Patches for which figures will be generated and saved as .eps files
selPatches = [4 5];
nSel = length(selPatches);


fName = fullfile(fiToolboxRootPath,'results','bootstrap','fl_Macbeth+Fl_bootstrap_100_1.mat');
flBoot = load(fName);
wave = flBoot.wave(:);

fName = fullfile(fiToolboxRootPath,'results','bootstrap','em_Macbeth+Fl_bootstrap_100_1.mat');
emBoot = load(fName);


sr = 1:5:length(wave);

%% Display data for all patches

% Pixel values

flMeasValsEst = cell2mat(shiftdim(flBoot.reflValsEst,-3)) + cell2mat(shiftdim(flBoot.flValsEst,-3));
avgFlMeasValsEst = mean(flMeasValsEst,4);
avgFlMeasVals = mean(flBoot.measVals,4);


figure;
hold all; grid on; box on;
plot(avgFlMeasValsEst(:),avgFlMeasVals(:),'.');
xlabel('Model predicted pixel value');
ylabel('Measured pixel value');


figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    
    
    tmp1 = avgFlMeasValsEst(:,:,sampleID)';
    tmp2 = avgFlMeasVals(:,:,sampleID)';
    
    plot(tmp1,tmp2,'.');

end
end


% Reflectance

flReflEst = cell2mat(shiftdim(flBoot.reflEst,-2));
avgFlReflEst = mean(flReflEst,3);
sorted = sort(flReflEst,3);
avgFlReflEstLb = sorted(:,:,3);
avgFlReflEstUb = sorted(:,:,97);


emReflEst = cell2mat(shiftdim(emBoot.reflEst,-2));
avgEmReflEst = mean(emReflEst,3);
sorted = sort(emReflEst,3);
avgEmReflEstLb = sorted(:,:,3);
avgEmReflEstUb = sorted(:,:,97);


figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot(wave,avgFlReflEst(:,sampleID),'g','LineWidth',2);
    plot(wave,avgEmReflEst(:,sampleID),'k','LineWidth',2);
    
    
    % Confidence intervals
    plot(wave,avgFlReflEstLb(:,sampleID),'c','LineWidth',1);
    plot(wave,avgFlReflEstUb(:,sampleID),'c','LineWidth',1);
    
    % Confidence intervals
    plot(wave,avgEmReflEstLb(:,sampleID),'k','LineWidth',1);
    plot(wave,avgEmReflEstUb(:,sampleID),'k','LineWidth',1);
    
    % Reference
    plot(wave,flBoot.reflRef(:,sampleID),'b--','LineWidth',2);
    
    
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);

    rmse = sqrt(mean((avgFlReflEst(:,sampleID) - flBoot.reflRef(:,sampleID)).^2));
    title(sprintf('RMSE %.2f',rmse));

end
end

% Emission (normalized)

flEmEst = cell2mat(shiftdim(flBoot.emEst,-2));
avgFlEmEst = mean(flEmEst,3);
sorted = sort(flEmEst,3);
avgFlEmEstLb = sorted(:,:,3);
avgFlEmEstUb = sorted(:,:,97);

emEst = cell2mat(shiftdim(emBoot.emEst,-2));
avgEmEst = mean(emEst,3);
sorted = sort(emEst,3);
avgEmEstLb = sorted(:,:,3);
avgEmEstUb = sorted(:,:,97);


figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    
    flNF = max(avgFlEmEst(:,sampleID));
    nF = max(avgEmEst(:,sampleID));
    
    plot(wave,avgFlEmEst(:,sampleID)/flNF,'m','LineWidth',2);
    
    plot(wave,avgFlEmEstLb(:,sampleID)/flNF,'c','LineWidth',1);
    plot(wave,avgFlEmEstUb(:,sampleID)/flNF,'c','LineWidth',1);
    
    
    plot(wave,avgEmEst(:,sampleID)/nF,'g','LineWidth',2);
    
    plot(wave,avgEmEstLb(:,sampleID)/nF,'k','LineWidth',1);
    plot(wave,avgEmEstUb(:,sampleID)/nF,'k','LineWidth',1);
    
    
    plot(wave,flBoot.emRef(:,sampleID)/max(flBoot.emRef(:,sampleID)),'b--','LineWidth',2);
    xlim([min(wave) max(wave)]);
    ylim([0.0 1.05]);

end
end

% Excitation

flExEst = cell2mat(shiftdim(flBoot.exEst,-2));


avgFlExEst = mean(flExEst,3);
sorted = sort(flExEst,3);
avgFlExEstLb = sorted(:,:,3);
avgFlExEstUb = sorted(:,:,97);


maxVal = max([avgFlEmEst(:); flBoot.emRef(:)]);
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot(wave,avgFlExEst(:,sampleID),'g','LineWidth',2);
    
    plot(wave,avgFlExEstLb(:,sampleID),'c','LineWidth',1);
    plot(wave,avgFlExEstUb(:,sampleID),'c','LineWidth',1);
    
    plot(wave,flBoot.exRef(:,sampleID),'b--','LineWidth',2);
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);

    rmse = sqrt(mean((avgFlExEst(:,sampleID) - flBoot.exRef(:,sampleID)).^2));
    title(sprintf('RMSE %.2f',rmse));

end
end



%% Now generate figures for specific, selected patches

% Pixel values

flMeasValsEst = cell2mat(shiftdim(flBoot.reflValsEst,-3)) + cell2mat(shiftdim(flBoot.flValsEst,-3));
avgFlMeasValsEst = mean(flMeasValsEst,4);
sorted = sort(flMeasValsEst,4);
avgFlMeasValsEstLb = sorted(:,:,:,3);
avgFlMeasValsEstUb = sorted(:,:,:,97);

emMeasValsEst = cell2mat(shiftdim(emBoot.reflValsEst,-3)) + cell2mat(shiftdim(emBoot.flValsEst,-3));
avgEmMeasValsEst = mean(emMeasValsEst,4);
sorted = sort(emMeasValsEst,4);
avgEmMeasValsEstLb = sorted(:,:,:,3);
avgEmMeasValsEstUb = sorted(:,:,:,97);


for s=1:nSel
   
    id = selPatches(s);
    tmp = flBoot.measVals(:,:,id); tmp = tmp(:);


    figure;
    hold on; grid on; box on;    
     
    % Single fluorophore estimate
    tmp2 = avgFlMeasValsEst(:,:,id);
    l(2) = plot(tmp(:),tmp2(:),'go','markerSize',ms,'lineWidth',lw);
    
    % CIM estimate
    tmp2 = avgEmMeasValsEst(:,:,id);
    l(3) = plot(tmp(:),tmp2(:),'kx','markerSize',ms,'lineWidth',lw);
    
    % Reference
    plot(linspace(0,1.1,10),linspace(0,1.1,10),'r-','lineWidth',0.5*lw);
    
    xlim([0 1.05]);
    ylim([0 1.05]);
    
    xlabel('Measured pixel intensity','fontsize',fs);
    ylabel('Predicted pixel intensity','fontsize',fs);
    
    legend(l(2:3),leg(2:3),'location','northwest');
    
    set(gca,'yTick',0:0.2:1);
    set(gca,'fontsize',fs-2);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPosition',sz);
    
    if ~isempty(saveDir)
        fName = fullfile(saveDir,sprintf('flPixels_%i.eps',id));
        print('-depsc',fName);
    end
end



%  Reflectance
for s=1:nSel
   
    id = selPatches(s);
    
    figure;
    hold on; grid on; box on;
    
    % Confidence interval on the CIM estimate
    fill([wave; rot90(wave,2)],[avgEmReflEstUb(:,id); rot90(avgEmReflEstLb(:,id),2)],[0.8 0.8 0.8],'edgeColor','none');
   
    
    % Confidence interval on single fluorophore estimate
    fill([wave; rot90(wave,2)],[avgFlReflEstUb(:,id); rot90(avgFlReflEstLb(:,id),2)],[0.8 1 0.8],'edgeColor','none');
    
     
    % Single fluorophore estimate
    l(2) = plot(wave,avgFlReflEst(:,id),'g--','lineWidth',lw);
    plot(wave(sr),avgFlReflEst(sr,id),'go','markerSize',ms,'lineWidth',lw*0.5);
    
    % CIM estimate
    l(3) = plot(wave,avgEmReflEst(:,id),'k--','lineWidth',lw);
    plot(wave(sr),avgEmReflEst(sr,id),'kx','markerSize',ms,'lineWidth',lw*0.5);
    
    % Reference
    l(1) = plot(wave,flBoot.reflRef(:,id),'r','lineWidth',lw);
    
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);
    
    lh = legend(l,leg,'location','northwest');
    
    ch = get(lh,'Children');
    set(ch(4),'Marker','o','MarkerSize',ms);
    set(ch(1),'Marker','x','MarkerSize',ms);
    
    xlabel('Wavelength, nm','fontsize',fs);
    set(gca,'yTick',0:0.2:1);
    set(gca,'fontsize',fs-2);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPosition',sz);
    
    if ~isempty(saveDir)
        fName = fullfile(saveDir,sprintf('flRefl_%i.eps',id));
        print('-depsc',fName);
    end
end


%  Emission (normalized)
for s=1:nSel
   
    id = selPatches(s);
    nF = max(avgFlEmEst(:,id));
    nFem = max(avgEmEst(:,id));
    
    figure;
    hold on; grid on; box on;
    
    % Confidence interval on CIM estimate
    fill([wave; rot90(wave,2)],[avgEmEstUb(:,id)/nFem; rot90(avgEmEstLb(:,id)/nFem,2)],[0.8 0.8 0.8],'edgeColor','none');
    
    % Confidence interval on single fluorophore estimate
    fill([wave; rot90(wave,2)],[avgFlEmEstUb(:,id)/nF; rot90(avgFlEmEstLb(:,id)/nF,2)],[0.8 1 0.8],'edgeColor','none');
    
    % Single fluorophore estimate
    l(2) = plot(wave,avgFlEmEst(:,id)/nF,'g--','lineWidth',lw);
    plot(wave(sr),avgFlEmEst(sr,id)/nF,'go','markerSize',ms,'lineWidth',lw*0.5);

    
    % Single emission estimate
    l(3) = plot(wave,avgEmEst(:,id)/nFem,'k--','lineWidth',lw);
    plot(wave(sr),avgEmEst(sr,id)/nFem,'kx','markerSize',ms,'lineWidth',lw*0.5);
    
    % Reference
    l(1) = plot(wave,flBoot.emRef(:,id)/max(flBoot.emRef(:,id)),'r','lineWidth',lw);
    
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);
    
    lh = legend(l,leg,'location','northeast');
    ch = get(lh,'Children');
    set(ch(4),'Marker','o','MarkerSize',ms);
    set(ch(1),'Marker','x','MarkerSize',ms);
    
    xlabel('Wavelength, nm','fontsize',fs);
    set(gca,'yTick',0:0.2:1);
    set(gca,'fontsize',fs-2);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPosition',sz);
    
    if ~isempty(saveDir)
        fName = fullfile(saveDir,sprintf('flEm_%i.eps',id));
        print('-depsc',fName);
    end
end

%  Excitation (normalized)
l = [];
for s=1:nSel
   
    id = selPatches(s);
    nF = max(avgFlExEst(:,id));
    
    figure;
    hold on; grid on; box on;
    
    % Confidence interval on single fluorophore estimate
    fill([wave; rot90(wave,2)],[avgFlExEstUb(:,id)/nF; rot90(avgFlExEstLb(:,id)/nF,2)],[0.8 1 0.8],'edgeColor','none');
    
    % Single fluorophore estimate
    l(2) = plot(wave,avgFlExEst(:,id)/nF,'g--','lineWidth',lw);
    plot(wave(sr),avgFlExEst(sr,id)/nF,'go','markerSize',ms,'lineWidth',lw*0.5);

    
    % Reference
    l(1) = plot(wave,flBoot.exRef(:,id)/max(flBoot.exRef(:,id)),'r','lineWidth',lw);
    
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);
    
    lh = legend(l(1:2),leg(1:2),'location','northeast');
    ch = get(lh,'Children');
    set(ch(1),'Marker','o','MarkerSize',ms);
    
    xlabel('Wavelength, nm','fontsize',fs);
    set(gca,'yTick',0:0.2:1);
    set(gca,'fontsize',fs-2);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPosition',sz);
    
    if ~isempty(saveDir)
        fName = fullfile(saveDir,sprintf('flEx_%i.eps',id));
        print('-depsc',fName);
    end
end


% Excitation (absolute)
l = [];
for s=1:nSel
   
    id = selPatches(s);
    nF = max(avgFlEmEst(:,id));
    
    figure;
    hold on; grid on; box on;
    
    % Confidence interval on single fluorophore estimate
    fill([wave; rot90(wave,2)],[avgFlExEstUb(:,id)*nF; rot90(avgFlExEstLb(:,id)*nF,2)],[0.8 1 0.8],'edgeColor','none');
    
    % Single fluorophore estimate
    l(2) = plot(wave,avgFlExEst(:,id)*nF,'g--','lineWidth',lw);
    plot(wave(sr),avgFlExEst(sr,id)*nF,'go','markerSize',ms,'lineWidth',lw*0.5);

    
    % Reference
    l(1) = plot(wave,flBoot.exRef(:,id)*max(flBoot.emRef(:,id)),'r','lineWidth',lw);
    
    
    xlim([min(wave) max(wave)]);
    
    lh = legend(l(1:2),leg(1:2),'location','northeast');
    ch = get(lh,'Children');
    set(ch(1),'Marker','o','MarkerSize',ms);
    
    xlabel('Wavelength, nm','fontsize',fs);
    set(gca,'fontsize',fs-2);
    set(gcf,'Units','centimeters');
    set(gcf,'PaperPosition',sz);
    
    if ~isempty(saveDir)
        fName = fullfile(saveDir,sprintf('flExToScale_%i.eps',id));
        print('-depsc',fName);
    end
end

