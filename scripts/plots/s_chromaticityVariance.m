close all;
clear all;
clc;

wave = 380:4:1000;

fName = fullfile(fiToolboxRootPath,'data','redFl');
redFl = fiReadFluorophore(fName,'wave',wave);

fName = fullfile(fiToolboxRootPath,'data','greenFl');
greenFl = fiReadFluorophore(fName,'wave',wave);

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');

fs = 10;
lw = 1;
ms = 10;

lw2 = 0.5;
sz = [1 1 5 1.5];
sz2 = [1 1 3.5 3.5];

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

set(gcf,'Units','centimeters');
set(gcf,'PaperPosition',sz);

xlabel('Wavelength, nm','fontsize',fs);
set(gca,'fontsize',fs-2)
lg = legend(l,{'$e_x(\lambda)$','$e_m(\lambda)$'},'fontsize',fs,'location','northeast','interpreter','latex');


% Change everything to black
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


fName = fullfile(saveDir,'TwoFluorophores.eps');
print('-depsc',fName);


%% Emission spectra
spd = zeros(length(wave));
for i=1:length(wave)
    
    light = zeros(length(wave),1);
    light(i) = 1;
    
    spd(:,i) = redFl.emission*redFl.excitation'*light + greenFl.emission*greenFl.excitation'*light;
    
end
spd = spd*diag(1./max(spd));


selWaveIDs = [10 40 60];

figure; 
hold on; grid off; box on;
plot(wave,spd(:,selWaveIDs),'LineWidth',lw);
xlim([min(wave) 800]);
ylim([-0.05 1.05]);

set(gca,'XTick',400:200:1000);
set(gca,'YTick',0:0.2:1);

set(gcf,'Units','centimeters');
set(gcf,'PaperPosition',sz);

set(gca,'fontsize',fs-2)
xlabel('Wavelength, nm','fontsize',fs);

for i=1:3
    legLabels{i} = sprintf('%inm',wave(selWaveIDs(i)));
end

legend(legLabels,'FontSize',fs,'location','northeast');




fName = fullfile(saveDir,'chrInv.eps');
print('-depsc',fName);




%% Chromaticity diagram.

% Compute spectra and convert to xyY
spd(isnan(spd)) = 0;
spd = spd*1e17;

XYZ = ieXYZFromPhotons(spd',wave);

x = XYZ(:,1)./sum(XYZ,2);
y = XYZ(:,2)./sum(XYZ,2);


% Plot on the CIE diagram
figure; 
cieplot();
plot(x,y,'kx','LineWidth',lw,'MarkerSize',ms);
xlabel('x','fontsize',fs);
ylabel('y','fontsize',fs);


set(gca,'XTick',0:0.2:0.8);
set(gca,'YTick',0:0.2:0.8);

set(gcf,'Units','centimeters');
set(gcf,'PaperPosition',sz2);
set(gca,'FontSize',fs-2);

fName = fullfile(saveDir,'chrxy.eps');
print('-depsc','-opengl','-r400',fName);




    