close all;
clear all;
clc;

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');


%% Single fluorophore algorithm
fName = fullfile(fiToolboxRootPath,'results','evaluation','McNamara-Boswell_simNBands_fl');
load(fName);

fs = 8;
sz = [1 1 8 8];
lw = 2;
cntrs = [0.1 0.075 0.05 0.025];

figure;
grid on; hold on;
hndl = surf(unique(filterGrid),unique(channelGrid),dMatErr);
xlabel('# camera filters','fontsize',fs);
ylabel('# illuminants','fontsize',fs);
zlabel('RMSE','fontsize',fs);
set(gca,'fontsize',fs-2);
set(gca,'ZScale','log');
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
set(gca,'ZMinorGrid','off');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
zlim([0.01, 0.5]);
xlim([min(filterGrid(:)) max(channelGrid(:))]);
ylim([min(filterGrid(:)) max(channelGrid(:))]);
campos([  169.8802  196.8408    1.6022]);

cMap = colormap(gca);

C = contourc(unique(filterGrid),unique(channelGrid),dMatErr,cntrs);
offset = 1;
for i=1:length(cntrs)
    lvl = C(1,offset);
    nPts = C(2,offset);
    clr = round(interp1(linspace(min(dMatErr(:)),max(dMatErr(:)),size(cMap,1)),1:size(cMap,1),lvl));
    
    plot3(C(1,offset+1:offset+nPts),C(2,offset+1:offset+nPts),ones(1,nPts)*0.01,'lineWidth',lw,'color',cMap(clr,:));
    offset = offset + nPts + 1;
end
uistack(hndl,'top');

fName = fullfile(saveDir,'flNBands.eps');
print('-depsc',fName);

%% Multi fluorophore algorithm

