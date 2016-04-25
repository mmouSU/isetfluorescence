close all;
clear all;
clc;

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');


fName = fullfile(fiToolboxRootPath,'results','evaluation','McNamara-Boswell_multiFl_nBasis');
load(fName);

fs = 8;
sz = [1 1 8 8];
lw = 2;
cntrs = [0.05 0.025 0.01 0.0075 0.005];

figure;
grid on; hold on;
hndl = surf(unique(emBasisGrid),unique(exBasisGrid),dMatErr);
xlabel('# emission basis','fontsize',fs);
ylabel('# excitation basis','fontsize',fs);
zlabel('RMSE','fontsize',fs);
set(gca,'fontsize',fs-2);
set(gca,'ZScale','log');
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
set(gca,'ZMinorGrid','off');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
zlim([0.003,0.2]);
xlim([min(emBasisGrid(:)) max(emBasisGrid(:))]);
ylim([min(exBasisGrid(:)) max(exBasisGrid(:))]);
campos([150.7973  174.4253    0.4096]);

cMap = colormap(gca);

C = contourc(unique(emBasisGrid),unique(exBasisGrid),dMatErr,cntrs);
offset = 1;
for i=1:length(cntrs)
    lvl = C(1,offset);
    nPts = C(2,offset);
    clr = round(interp1(linspace(min(dMatErr(:)),max(dMatErr(:)),size(cMap,1)),1:size(cMap,1),lvl));
    
    plot3(C(1,offset+1:offset+nPts),C(2,offset+1:offset+nPts),ones(1,nPts)*0.003,'lineWidth',lw,'color',cMap(clr,:));
    offset = offset + nPts + 1;
end
uistack(hndl,'top');

fName = fullfile(saveDir,'multiFlNBasis.eps');
print('-depsc',fName);