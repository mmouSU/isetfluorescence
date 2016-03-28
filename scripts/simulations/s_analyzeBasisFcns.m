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
lw = 2;
mkSz = 5;
fs = 6;
figSize = [0 0 2.25 1.5];

figure;
hold all; grid on; box on;
pl(1) = plot(1:24,cumNormReflScore,lineStyle{2});
pl(2) = plot(1:24,cumNormExScore,lineStyle{3});
pl(3) = plot(1:24,cumNormEmScore,lineStyle{4});
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
legend({'Reflectance','Excitation','Emission'},'location','southeast');
ylim([0.5 1.1]);
xlabel('# of basis functions');
ylabel('Variance explained');
set(gca,'fontSize',fs);
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
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
xlabel('Wavelength, nm');
xlim([min(wave) max(wave)]);
set(gca,'fontSize',fs);
print('-depsc',fullfile(saveDir,'reflBasis.eps'));

% Excitation basis

figure;
hold all; grid on; box on;
pl = plot(wave,exBasis(:,1:3));
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
xlabel('Wavelength, nm');
xlim([min(wave) max(wave)]);
set(gca,'fontSize',fs);
print('-depsc',fullfile(saveDir,'exBasis.eps'));

% Emission basis

figure;
hold all; grid on; box on;
pl = plot(wave,emBasis(:,1:3));
set(pl,'lineWidth',lw,'markerSize',mkSz);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
xlabel('Wavelength, nm');
xlim([min(wave) max(wave)]);
set(gca,'fontSize',fs);
print('-depsc',fullfile(saveDir,'emBasis.eps'));




