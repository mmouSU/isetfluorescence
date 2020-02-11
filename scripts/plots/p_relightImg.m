% Render images of a test target under different illuminants. This script
% generates synthetic images shown in the main manuscript as well as
% the Supplemental Material.
%
% Note, run p_relightImg_Natural_Comp.m to run a comparison with real captured
% data.
%
% Copyright, Henryk Blasinski

close all;
clear all;
clc;

ieInit;
wave = 380:4:1000;
nWaves = length(wave);

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
saveDir = fullfile('~','Desktop','Test images');
% saveDir = [];
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

dataDir = fullfile(fiToolboxRootPath,'results','experiments');

targetName = 'Natural';
% targetName = 'Target6';

switch targetName
    case 'Natural'
        nCols = 128;
        nRows = 103;
        
        % Load illuminants
        fName = fullfile(fiToolboxRootPath,'camera','illuminants');
        illuminant = ieReadSpectra(fName,wave);
        illSubset = illuminant(:,[2 6]);

        fName = fullfile(fiToolboxRootPath,'data','Broadband9800K');
        illSubset = [illSubset, ieReadSpectra(fName,wave)]; 
        fName = fullfile(fiToolboxRootPath,'data','Broadband6500K');
        illSubset = [illSubset, ieReadSpectra(fName,wave)];
        fName = fullfile(fiToolboxRootPath,'data','Broadband2000K');
        illSubset = [illSubset, ieReadSpectra(fName,wave)];
        fName = fullfile(fiToolboxRootPath,'data','NarrowbandBlue');
        illSubset = [illSubset, ieReadSpectra(fName,wave)];
        

    case 'Target6'
        nCols = 128;
        nRows = 103;
        
        % Load illuminants
        fName = fullfile(fiToolboxRootPath,'data',targetName,'LEDCubeRadiance','illuminants.mat');
        ledCubeLights = ieReadSpectra(fName,wave);
        illSubset = ledCubeLights(:,10);
        
        fName = fullfile(fiToolboxRootPath,'data','Broadband9800K');
        illSubset = [illSubset, ieReadSpectra(fName,wave)]; 
        
    otherwise
        error('Define illuminants for target %s',targetName);
end

 

% Load scene reflectance and fluorescence
reflArray = zeros(nRows,nCols,156);
for cc=1:nCols
    
    fName = fullfile(dataDir,sprintf('multiFl_%s_col_%i.mat',targetName,cc));
    data = load(fName);
    
    reflArray(:,cc,:) = data.reflEst';
    
    for rr=1:nRows
        deltaL = data.wave(2) - data.wave(1);
        dm = data.dMatEst{rr}/deltaL;
        dm(wave >= 650,wave >= 650) = 0;
        flArray(rr,cc) = fluorophoreCreate('type','fromdonaldsonmatrix',...
                                         'DonaldsonMatrix',dm,...
                                         'wave',data.wave);
    end
    
end    

%% Create a scene with the specified reflectance properties and a particular
% illuminant
sceneReTemplate = sceneReflectanceArray(reflArray,1,wave);
sceneFlTemplate = fluorescentSceneCreate('type','fromfluorophore','fluorophore',flArray);

fName = fullfile(fiToolboxRootPath,'data','macbethChart');
refl = ieReadSpectra(fName,wave);
sceneMacbeth = sceneReflectanceArray(reshape(refl',[4 6 nWaves]),10,wave);

% Camera filter transmissivities
fName = fullfile(fiToolboxRootPath,'data','CanonG7x');
cfa = ieReadSpectra(fName,wave);


%%

for i=1:size(illSubset,2)      
        
    sceneMacbeth = sceneAdjustIlluminant(sceneMacbeth,illSubset(:,i),0);
    ieAddObject(sceneMacbeth);
    
    sceneRe = sceneAdjustIlluminant(sceneReTemplate,illSubset(:,i),0);
    sceneRe = sceneSet(sceneRe,'name',sprintf('refl - ill%i',i));
    ieAddObject(sceneRe);
    
    sceneFl = fiSceneAddFluorescence(sceneRe,sceneFlTemplate,'replace',true);
    sceneFl = sceneSet(sceneFl,'name',sprintf('fl - ill%i',i));
    ieAddObject(sceneFl);
    
    % Add fluorescence
    sceneReFl = fiSceneAddFluorescence(sceneRe,sceneFlTemplate);
    sceneReFl = sceneSet(sceneReFl,'name',sprintf('refl+fl - ill%i',i));
    ieAddObject(sceneReFl);
    
    
    sceneWindow;
    
    % Simulate the entire camera pipeline
    oi = oiCreate;
    oiRe = oiCompute(oi,sceneRe);
    oiFl = oiCompute(oi,sceneFl);
    oiReFl = oiCompute(oi,sceneReFl);
    oiMacbeth = oiCompute(oi,sceneMacbeth);
    
    sensor = sensorCreate;
    sensor = sensorSet(sensor,'wave',wave);
    sensor = sensorSet(sensor,'filter transmissivities',cfa);
    sensor = sensorSetSizeToFOV(sensor,[sceneGet(sceneRe,'fov horizontal') sceneGet(sceneRe,'fov vertical')],sceneRe,oi);
    
    sensorMacbeth = sensorCompute(sensor,oiMacbeth);

    
    sensorReFl = sensorCompute(sensor,oiReFl);
    expTime = sensorGet(sensorReFl,'exposure time');
    
    % We want all the images to be exposed the same way. The
    % reflectance+fluorescence image will always have the strongest signal.
    sensor = sensorSet(sensor,'exposure time',expTime);
    sensor = sensorSet(sensor,'auto exposure',0);
    
    sensorRe = sensorCompute(sensor,oiRe);
    sensorFl = sensorCompute(sensor,oiFl);
    
    ip = ipCreate;
    ipRe = ipCompute(ip,sensorRe);
    ipRe = ipSet(ipRe,'name',sprintf('refl - ill%i',i));
    
    ipFl = ipCompute(ip,sensorFl);
    ipFl = ipSet(ipFl,'name',sprintf('fl - ill%i',i));
    ipReFl = ipCompute(ip,sensorReFl);
    ipReFl = ipSet(ipReFl,'name',sprintf('refl+fl - ill%i',i));

    ipMacbeth = ipCompute(ip,sensorMacbeth);
    
    lRGBRe = ipGet(ipRe,'sensor channels');
    lRGBFl = ipGet(ipFl,'sensor channels');
    lRGBReFl = ipGet(ipReFl,'sensor channels');
    lRGBMacbeth = ipGet(ipMacbeth,'sensor channels');
    
    sRGBRe = ipGet(ipRe,'data srgb');
    sRGBFl = ipGet(ipFl,'data srgb');
    sRGBReFl = ipGet(ipReFl,'data srgb');
    sRGBMacbeth = ipGet(ipMacbeth,'data srgb');
    
    % Save images
    if ~isempty(saveDir)
        fName = fullfile(saveDir,sprintf('%s_renderedLight_%i_re.png',lower(targetName), i));
        imwrite(lRGBRe,fName);
        fName = fullfile(saveDir,sprintf('%s_renderedLight_%i_fl.png',lower(targetName), i));
        imwrite(lRGBFl,fName);
        fName = fullfile(saveDir,sprintf('%s_renderedLight_%i_reFl.png',lower(targetName), i));
        imwrite(lRGBReFl,fName);
        fName = fullfile(saveDir,sprintf('%s_renderedLight_%i_Macbeth.png',lower(targetName), i));
        imwrite(lRGBMacbeth,fName);
    end
    
    linImg = [lRGBRe lRGBFl lRGBReFl lRGBMacbeth];
    linImg = linImg./max(linImg(:));
    
    figure; imshow(linImg.^(1/2.2));
    
    ieAddObject(ipRe);
    ieAddObject(ipFl);
    ieAddObject(ipReFl);

    ipWindow;
    drawnow;
    
end