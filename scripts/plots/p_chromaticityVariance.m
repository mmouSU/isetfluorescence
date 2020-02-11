% This script demonstrates the variance of fluorescence emission when a
% mixture of fluorophores is present and generates Fig. 3 in the paper.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
% saveDir = fullfile('~','Desktop','Figures');
saveDir = [];
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

wave = 380:4:1000;

% Figure style parameters
fs = 10;
lw = 1;
ms = 10;

lw2 = 0.5;
sz = [1 1 5 1.5]*2.5;
sz2 = [1 1 3.5 3.5]*2.5;

markers = {'o','x','^'};
colors = {'r','g','b'};
selWaves = 1:3:length(wave);




% Load the fluorophores used in experiments
fName = fullfile(fiToolboxRootPath,'data','redFl');
redFl = fiReadFluorophore(fName,'wave',wave);

fName = fullfile(fiToolboxRootPath,'data','greenFl');
greenFl = fiReadFluorophore(fName,'wave',wave);



%% Plot excitation and emission spectra
figure;
hold on; grid off; box on;

l(1) = plot(wave,redFl.excitation,'r--','lineWidth',lw);
plot(wave,greenFl.excitation,'g--','LineWidth',lw);
l(2) = plot(wave,redFl.emission/max(redFl.emission),'r-','lineWidth',lw);
plot(wave,greenFl.emission/max(greenFl.emission),'g-','LineWidth',lw);

xlim([min(wave) 800]);
ylim([-0.05 1.05]);
set(gca,'XTick',400:200:1000);
set(gca,'YTick',0:0.2:1);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'Units','centimeters');
set(gcf,'PaperPosition',sz);
xlabel('Wavelength, nm','fontsize',fs,'interpreter','latex');
set(gca,'fontsize',fs-2)
lg = legend(l,{'$e_x(\lambda)$','$e_m(\lambda)$'},'fontsize',fs,'location','northeast','interpreter','latex');


% Change all legend symbols to black
ch = get(lg,'Children');
set(ch,'Color','black');

lineHndl = findobj(ch,'Type','line');
line = get(lineHndl,'XData');
for i=1:length(lineHndl)
    if length(line{i}) == 2
        line{i}(2) = line{i}(2)*3/4;
        set(lineHndl(i),'XData',line{i});
    end
end

textHndl = findobj(ch,'Type','text');
textPos = get(textHndl,'Position');
for i=1:length(textHndl)
    textPos{i}(1) = textPos{1}(1);
    set(textHndl(i),'Position',textPos{i});
end

if ~isempty(saveDir)
    fName = fullfile(saveDir,'twoFluorophores.eps');
    print('-depsc',fName);
end


%% Plot the (normalized) emission radiances of the two fluorophore composition

% We assume that the sample is illuminated with monochromatic light and we
% select radiances at three different excitation wavelengths.
selWaveIDs = [10 40 60];

spd = zeros(length(wave));
for i=1:length(wave)
    
    light = zeros(length(wave),1);
    light(i) = 1;
    
    spd(:,i) = redFl.emission*redFl.excitation'*light + greenFl.emission*greenFl.excitation'*light;
    
end
spd = spd*diag(1./max(spd));



figure; 
hold all; grid off; box on;

for i=1:length(selWaveIDs)
    plot(wave,spd(:,selWaveIDs(i)),colors{i},'LineWidth',lw);
    lg(i) = plot(wave(selWaves),spd(selWaves,selWaveIDs(i)),[markers{i} colors{i}],'LineWidth',lw,'markersize',ms/2);
end
xlim([min(wave) 800]);
ylim([-0.05 1.05]);

set(gca,'XTick',400:200:1000);
set(gca,'YTick',0:0.2:1);
set(gca,'TickLabelInterpreter','latex');
set(gcf,'Units','centimeters');
set(gcf,'PaperPosition',sz);

set(gca,'fontsize',fs-2)
xlabel('Wavelength, nm','fontsize',fs,'interpreter','latex');

for i=1:3
    legLabels{i} = sprintf('%inm',wave(selWaveIDs(i)));
end

lh = legend(lg,legLabels,'FontSize',fs,'location','northeast','interpreter','latex');

if ~isempty(saveDir)
    fName = fullfile(saveDir,'chrInv.eps');
    print('-depsc',fName);
end



%% Chromaticity diagram.

% Compute spectra and convert to xyY. Multiply unit intensity spectra by a
% large number so that they represent plausible photon units.
spd(isnan(spd)) = 0;
spd = spd*1e17;

XYZ = ieXYZFromPhotons(spd',wave);

x = XYZ(:,1)./sum(XYZ,2);
y = XYZ(:,2)./sum(XYZ,2);


% Plot on the CIE diagram
figure; 
cieplot();
set(gca,'FontSize',fs-2);
plot(x,y,'kx','LineWidth',lw,'MarkerSize',ms);
xlabel('x','fontsize',fs,'interpreter','latex');
ylabel('y','fontsize',fs,'interpreter','latex');

set(gca,'TickLabelInterpreter','latex');
set(gca,'XTick',0:0.2:0.8);
set(gca,'YTick',0:0.2:0.8);

set(gcf,'Units','centimeters');
set(gcf,'PaperPosition',sz2);


if ~isempty(saveDir)
    fName = fullfile(saveDir,'chrxy.eps');
    print('-depsc','-opengl','-r400',fName);
end



    