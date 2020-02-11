% Generate plots showing the accuracy of the multi- and single fluorophore
% estimation algorithms as a function of the SNR.
% This script reproduces Fig. 7 from the paper.
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

if ~exist(saveDir,'dir') && ~isempty(saveDir)
    mkdir(saveDir);
end

dataset = 'McNamara-Boswell';

lw = 1;
fs = 8;
ms = 3;
sz = [1 1 8 5];

%% Single fluorophore

fName = fullfile(fiToolboxRootPath,'results','evaluation',sprintf('%s_simSNR_Fl.mat',dataset));
load(fName);

figure;
hold on; grid on; box on;
errorbar(SNRdB,avgPixelErr,stdPixelErr,'-rs','lineWidth',lw,'markerSize',ms);
errorbar(SNRdB,avgReflErr,stdReflErr,'-gd','lineWidth',lw,'markerSize',ms);
errorbar(SNRdB,avgExErr,stdExErr,'-bo','lineWidth',lw,'markerSize',ms);
errorbar(SNRdB,avgEmErr,stdEmErr,'-c^','lineWidth',lw,'markerSize',ms);
% set(gca,'yscale','log');
set(gca,'fontsize',fs-2);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',sz);
set(gca,'XMinorGrid','off');
set(gca,'YMinorGrid','off');
set(gca,'TickLabelInterpreter','latex');
xlim([-15 max(SNRdB)]);
ylim([0 0.6]);
xlabel('SNR, dB','fontsize',fs,'Interpreter','LaTeX');
ylabel('RMSE','fontsize',fs,'Interpreter','LaTeX');
legend({'Pixel','Reflectance','Excitation','Emission'},'location','northeast','Interpreter','LaTeX');

if ~isempty(saveDir)
    fName = fullfile(saveDir,'flSNR.eps');
    print('-depsc',fName);
end

%% Multi-fluorophore

fName = fullfile(fiToolboxRootPath,'results','evaluation',sprintf('%s_simSNR_multiFl.mat',dataset));
load(fName);

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
set(gca,'TickLabelInterpreter','latex');
xlim([-15 max(SNRdB)]);
ylim([0 0.6]);
xlabel('SNR, dB','fontsize',fs,'Interpreter','LaTeX');
ylabel('RMSE','fontsize',fs,'Interpreter','LaTeX');
legend({'Pixel','Reflectance','Donaldson matrix'},'location','northeast','Interpreter','LaTeX');

if ~isempty(saveDir)
    fName = fullfile(saveDir,'multiFlSNR.eps');
    print('-depsc',fName);
end



