close all;
clear all;
clc;

wave = 380:4:1068;

nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = createBasisSet('reflectance','wave',wave','n',24);
[exBasis, exScore] = createBasisSet('excitation','wave',wave','n',24);
[emBasis, emScore] = createBasisSet('emission','wave',wave','n',24);

cumNormReflScore = cumsum(reflScore/sum(reflScore));
cumNormExScore = cumsum(exScore/sum(exScore));
cumNormEmScore = cumsum(emScore/sum(emScore));

%% Generate plots
close all;

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');
lineStyle = {'rs-','gd-','bo-','c^-'};
lw = 1;
mkSz = 2;
fs = 5;
figSize = [1 1 4.5 3];

figure;
hold all; grid on; box on;
pl(1) = plot(1:24,cumNormReflScore,lineStyle{2});
pl(2) = plot(1:24,cumNormExScore,lineStyle{3});
pl(3) = plot(1:24,cumNormEmScore,lineStyle{4});
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',figSize);
legend({'Reflectance','Excitation','Emission'},'location','southeast');
ylim([0.5 1.1]);
xlabel('# of basis functions','fontsize',fs);
ylabel('Energy','fontsize',fs);
set(gca,'fontSize',fs-2);
print('-depsc',fullfile(saveDir,'basisVariance.eps'));

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
set(gcf,'PaperPosition',figSize);
set(gca,'XTick',[400 600 800 1000]);
set(gca,'YTick',-0.2:0.1:0.2);
xlabel('Wavelength, nm','fontsize',fs);
xlim([min(wave) max(wave)]);
ylim(1.05*[min(min(reflBasis(:,1:3))) max(max(reflBasis(:,1:3)))]);
set(gca,'fontSize',fs-2);
print('-depsc',fullfile(saveDir,'reflBasis.eps'));

% Excitation basis

figure;
hold all; grid on; box on;
pl = plot(wave,exBasis(:,1:3));
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',figSize);
xlabel('Wavelength, nm','fontsize',fs);
set(gca,'XTick',[400 600 800 1000]);
set(gca,'YTick',-0.2:0.1:0.2);
xlim([min(wave) max(wave)]);
ylim(1.05*[min(min(exBasis(:,1:3))) max(max(exBasis(:,1:3)))]);
set(gca,'fontSize',fs-2);
print('-depsc',fullfile(saveDir,'exBasis.eps'));

% Emission basis

figure;
hold all; grid on; box on;
pl = plot(wave,emBasis(:,1:3));
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',figSize);
xlabel('Wavelength, nm','fontsize',fs);
set(gca,'XTick',[400 600 800 1000]);
set(gca,'YTick',-0.2:0.1:0.2);
xlim([min(wave) max(wave)]);
ylim(1.05*[min(min(emBasis(:,1:3))) max(max(emBasis(:,1:3)))]);
set(gca,'fontSize',fs-2);
print('-depsc',fullfile(saveDir,'emBasis.eps'));




