% Render images of a test target under different illuminants. This script
% generates radiance plots from Fig. 10 in the manuscript.
%
% Copyright, Henryk Blasinski
%
close all;
clear all;
clc;

ieInit;
wave = 380:4:1000;
nWaves = length(wave);

%%
% Figure display parameters
fs = 10;
lw = 1;
ms = 4;
sz2 = [1 1 4 4] * 2;
sr = 1:5:length(wave);

leg = {'Reference','Estimate'};

targetName = 'Target6';

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
saveDir = fullfile('~','Desktop','TestFolder');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

%% Select illuminants

fName = fullfile(fiToolboxRootPath,'data',targetName,'LEDCubeRadiance','illuminants');
illSubset = ieReadSpectra(fName,wave);
ids = [10];
illSubset = illSubset(:,ids);

fName = fullfile(fiToolboxRootPath,'data','Broadband9800K');
illSubset = [illSubset, ieReadSpectra(fName,wave)]; 

illNames = {'NarrowbandBlue','Broadband9800K'};

%% Load ground truth data

% Load the test target reflectance
fName = fullfile(fiToolboxRootPath,'data','Target6','RedReflectance');
reflRef(:,1) = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'data','Target6','MixReflectance');
reflRef(:,2) = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'data','Target6','GreenReflectance');
reflRef(:,3) = ieReadSpectra(fName,wave);

reflArray = permute(reflRef,[3 2 1]);

% Load test target fluorescence
fName = fullfile(fiToolboxRootPath,'data','Target6','RedFluorophore');
redFl = fiReadFluorophore(fName,'wave',wave);

fName = fullfile(fiToolboxRootPath,'data','Target6','GreenFluorophore');
greenFl = fiReadFluorophore(fName,'wave',wave);

% Assign fluorophores to patches
fluorophores(1,1,1) = redFl;
fluorophores(1,1,2) = redFl;

fluorophores(1,2,1) = redFl;
fluorophores(1,2,2) = greenFl;

fluorophores(1,3,1) = greenFl;
fluorophores(1,3,2) = greenFl;

% Assign quantum efficiencies
flQe(1,1,1) = 0.5;
flQe(1,1,2) = 0.5;

flQe(1,2,1) = 0.4;
flQe(1,2,2) = 0.09;

flQe(1,3,1) = 0.5;
flQe(1,3,2) = 0.5;

% Create reference fluorescence only and reflectance only scenes
sceneFlReference = fluorescentSceneCreate('type','fromfluorophore','fluorophore',...
                                        fluorophores,'wave',wave,'qe',flQe);
sceneReReference = sceneReflectanceArray(reflArray,1,wave);


%% Load estimated data
dataDir = fullfile(fiToolboxRootPath,'results','experiments');
nCols = 128;
nRows = 103; 

% Load scene reflectance and fluorescence
reflArray = zeros(nRows,nCols,156);
for cc=1:nCols
    
    fName = fullfile(dataDir,sprintf('multiFl_%s_col_%i.mat',targetName,cc));
    data = load(fName);
    
    reflArray(:,cc,:) = data.reflEst';
    
    for rr=1:nRows
        deltaL = data.wave(2) - data.wave(1);
        dm = data.dMatEst{rr}/deltaL;
        flArray(rr,cc) = fluorophoreCreate('type','fromdonaldsonmatrix',...
                                         'DonaldsonMatrix',dm,...
                                         'wave',data.wave);
    end
    
end    


% Create estimated fluorescence only and reflectance only scenes
sceneReTemplate = sceneReflectanceArray(reflArray,1,wave);
sceneFlTemplate = fluorescentSceneCreate('type','fromfluorophore','fluorophore',flArray);

% Segment out patches in the estimated data.
patch(1).posEst = [41, 70];
patch(1).posRef = [1, 1];
patch(1).name = 'Red';

patch(2).posEst = [96, 89];
patch(2).posRef = [2, 1];
patch(2).name = 'Mix';

patch(3).posEst = [115, 44];
patch(3).posRef = [3, 1];
patch(3).name = 'Green';

delta = 10;

for i=1:size(illSubset,2)      
        
    sceneRe = sceneAdjustIlluminant(sceneReTemplate,illSubset(:,i),0);
    sceneRe.data.photons = zeros(size(sceneRe.data.photons));
    sceneRe = sceneSet(sceneRe,'name',sprintf('refl - ill%i',i));
    ieAddObject(sceneRe);
    
    % Add fluorescence
    sceneReFl = fiSceneAddFluorescence(sceneRe,sceneFlTemplate);
    sceneReFl = sceneSet(sceneReFl,'name',sprintf('refl+fl - ill%i',i));
    ieAddObject(sceneReFl);
    

    % Reference scene.
    
    sceneReRef = sceneAdjustIlluminant(sceneReReference,illSubset(:,i),0);
    sceneReRef.data.photons = zeros(size(sceneReRef.data.photons));
    sceneReRef = sceneSet(sceneReRef,'name',sprintf('refl - ill%i',i));
    ieAddObject(sceneReRef);
    
    % Add fluorescence
    sceneReFlRef = fiSceneAddFluorescence(sceneReRef,sceneFlReference);
    sceneReFlRef = sceneSet(sceneReFlRef,'name',sprintf('refl+fl - ill%i',i));
    ieAddObject(sceneReFlRef);
    
    sceneWindow;
    
    % Draw plots
    for p=1:length(patch)
        
        estEm = vcGetROIData(sceneReFl,[patch(p).posEst, delta, delta],'photons');
        estEm = mean(estEm,1);
        estEm = estEm(:)/max(estEm);
        
        
        measEm = vcGetROIData(sceneReFlRef,[patch(p).posRef 0 0],'photons');
        measEm = mean(measEm,1);
        measEm = measEm(:)/max(measEm);
        
        figure;
        hold on; grid on; box on;
        plot(wave,measEm,'r','lineWidth',lw);
        plot(wave,estEm,'g--','lineWidth',lw);
        plot(wave(sr),estEm(sr),'go','markerSize',ms,'lineWidth',lw*0.5);
        
        xlim([min(wave) max(wave)]);
        ylim([-0.05 1.05]);
        
        xlabel('Wavelength, nm','fontsize',fs, 'interpreter','latex');
        ylabel('Fluoresced radiance $\rho_f$','fontsize',fs, 'interpreter','latex');
        set(gca,'yTick',0:0.2:1);
        set(gca,'fontsize',fs-2);
        set(gca,'TickLabelInterpreter','latex');
        set(gcf,'Units','centimeters');
        set(gcf,'PaperPosition',sz2);
        
        [~, ch] = legend(leg,'location','northeast','interpreter','latex');
        set(ch(end), 'Marker','o','MarkerSize',ms);
        
        if ~isempty(saveDir)
            fName = fullfile(saveDir,sprintf('%s_radiance_%s_%s.eps',lower(targetName), patch(p).name, illNames{i}));
            print('-depsc',fName);
        end
        
    end
    
end




