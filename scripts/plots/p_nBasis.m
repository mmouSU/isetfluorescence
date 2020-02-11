% Generate plots showing the accuracy of the multi- and single fluorophore
% estimation algorithms as a function of the number of excitation and 
% emission basis functions. This script reproduces Fig. 5 from the paper.
%
% Copyright, Henryk Blasinski 2016.

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir =[], then
% figures are not saved.
% saveDir = fullfile('~','Desktop','Figures');
saveDir = [];
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

dataset = 'McNamara-Boswell';


fs = 8;
sz = [1 1 8 8];
lw = 2;
cntrs = [0.025 0.01 0.0075];

%% Single fluorophore

fName = fullfile(fiToolboxRootPath,'results','evaluation',sprintf('%s_simNBasis_fl',dataset));
load(fName);

figure;
grid on; hold on;
hndl = surf(unique(emBasisGrid),unique(exBasisGrid),dMatErr);
xlabel('\# emission basis','fontsize',fs,'interpreter','latex');
ylabel('\# excitation basis','fontsize',fs,'interpreter','latex');
zlabel('RMSE','fontsize',fs,'interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
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

if ~isempty(saveDir)
    fName = fullfile(saveDir,'flNBasis.eps');
    print('-depsc',fName);
end


%% Multi-fluorophore

fName = fullfile(fiToolboxRootPath,'results','evaluation',sprintf('%s_simNBasis_multiFl',dataset));
load(fName);

figure;
grid on; hold on;
hndl = surf(unique(emBasisGrid),unique(exBasisGrid),dMatErr);
xlabel('\# emission basis','fontsize',fs,'interpreter','latex');
ylabel('\# excitation basis','fontsize',fs,'interpreter','latex');
zlabel('RMSE','fontsize',fs,'interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
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

if ~isempty(saveDir)
    fName = fullfile(saveDir,'multiFlNBasis.eps');
    print('-depsc',fName);
end