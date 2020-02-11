% This script analyzes the accuracy of linear approximation of spectral
% quantities (reflectance, excitation and emission) and generates Fig. 5
% from the paper.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir =[], then
% figures are not saved.
saveDir = fullfile('~','Desktop','Figures');
% saveDir = [];
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

% Figure style parameters
lineStyle = {'rs-','gd-','bo-','c^-'};
lw = 1;
mkSz = 5;
fs = 10;
figSize = [1 1 4.5 3]*2.5;

wave = 380:4:1068;


% Generate linear basis functions for reflectance, excitation and emission
nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = fiCreateBasisSet('reflectance','wave',wave','n',24);
[exBasis, exScore] = fiCreateBasisSet('excitation','wave',wave','n',24);
[emBasis, emScore] = fiCreateBasisSet('emission','wave',wave','n',24);

cumNormReflScore = cumsum(reflScore/sum(reflScore));
cumNormExScore = cumsum(exScore/sum(exScore));
cumNormEmScore = cumsum(emScore/sum(emScore));

%% Generate plots

figure;
hold all; grid on; box on;
pl(1) = plot(1:24,cumNormReflScore,lineStyle{2});
pl(2) = plot(1:24,cumNormExScore,lineStyle{3});
pl(3) = plot(1:24,cumNormEmScore,lineStyle{4});
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',figSize);
set(gca,'TickLabelInterpreter','latex');
legend({'Reflectance','Excitation','Emission'},'location','southeast','interpreter','latex');
ylim([0.5 1.1]);
xlabel('\# of basis functions','fontsize',fs,'interpreter','latex');
ylabel('Energy','fontsize',fs,'interpreter','latex');
set(gca,'fontSize',fs-2);

if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'basisVariance.eps'));
end

critReflVar = interp1(1:24,cumNormReflScore,nReflBasis);
critExVar = interp1(1:24,cumNormExScore,nExBasis);
critEmVar = interp1(1:24,cumNormEmScore,nEmBasis);

fprintf('%i reflectance basis account for %f variance\n',nReflBasis,critReflVar);
fprintf('%i excitation basis account for %f variance\n',nExBasis,critExVar);
fprintf('%i emission basis account for %f variance\n',nEmBasis,critEmVar);


% Reflectance basis

figure;
hold all; grid on; box on;
pl = plot(wave,reflBasis(:,1:3));
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','centimeters');
set(gca,'TickLabelInterpreter','latex');
set(gcf,'PaperPosition',figSize);
set(gca,'XTick',[400 600 800 1000]);
set(gca,'YTick',-0.2:0.1:0.2);
xlabel('Wavelength, nm','fontsize',fs,'interpreter','latex');
xlim([min(wave) max(wave)]);
ylim(1.05*[min(min(reflBasis(:,1:3))) max(max(reflBasis(:,1:3)))]);
set(gca,'fontSize',fs-2);
if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'reflBasis.eps'));
end

% Excitation basis

figure;
hold all; grid on; box on;
pl = plot(wave,exBasis(:,1:3));
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',figSize);
set(gca,'TickLabelInterpreter','latex');
xlabel('Wavelength, nm','fontsize',fs,'interpreter','latex');
set(gca,'XTick',[400 600 800 1000]);
set(gca,'YTick',-0.2:0.1:0.2);
xlim([min(wave) max(wave)]);
ylim(1.05*[min(min(exBasis(:,1:3))) max(max(exBasis(:,1:3)))]);
set(gca,'fontSize',fs-2);
if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'exBasis.eps'));
end

% Emission basis

figure;
hold all; grid on; box on;
pl = plot(wave,emBasis(:,1:3));
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',figSize);
set(gca,'TickLabelInterpreter','latex');
xlabel('Wavelength, nm','fontsize',fs,'interpreter','latex');
set(gca,'XTick',[400 600 800 1000]);
set(gca,'YTick',-0.2:0.1:0.2);
xlim([min(wave) max(wave)]);
ylim(1.05*[min(min(emBasis(:,1:3))) max(max(emBasis(:,1:3)))]);
set(gca,'fontSize',fs-2);
if ~isempty(saveDir)
    print('-depsc',fullfile(saveDir,'emBasis.eps'));
end



