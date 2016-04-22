close all;
clear all;
clc;

ieInit;
wave = 380:4:1000;

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures','Canon G7x');

fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);

illSubset = illuminant(:,[1 2 3]);


tmp = illuminantCreate('blackbody',wave,10000);
illSubset = [illSubset illuminantGet(tmp,'energy')];
tmp = illuminantCreate('blackbody',wave,6500);
illSubset = [illSubset illuminantGet(tmp,'energy')];
tmp = illuminantCreate('blackbody',wave,2000);
illSubset = [illSubset illuminantGet(tmp,'energy')];


dataDir = fullfile(fiToolboxRootPath,'results','experiments');
nCols = 128;
nRows = 103; 

% Load scene reflectance
reflArray = zeros(nRows,nCols,156);
for cc=1:nCols
    
    fName = fullfile(dataDir,sprintf('multifl2_Natural_col_%i.mat',cc));
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

%%

for i=4:size(illSubset,2)      
        
    sceneRe = sceneAdjustIlluminant(sceneReTemplate,illSubset(:,i),0);
    sceneRe = sceneSet(sceneRe,'name',sprintf('refl - ill%i',i));
    vcAddObject(sceneRe);
    
    sceneFl = fiSceneAddFluorescence(sceneRe,sceneFlTemplate,'replace',true);
    sceneFl = sceneSet(sceneFl,'name',sprintf('fl - ill%i',i));
    vcAddObject(sceneFl);
    
    % Add fluorescence
    sceneReFl = fiSceneAddFluorescence(sceneRe,sceneFlTemplate);
    sceneReFl = sceneSet(sceneReFl,'name',sprintf('refl+fl - ill%i',i));
    vcAddObject(sceneReFl);
    
    
    sceneWindow;
    
    % Simulate the entire camera pipeline
    oi = oiCreate;
    oiRe = oiCompute(oi,sceneRe);
    oiFl = oiCompute(oi,sceneFl);
    oiReFl = oiCompute(oi,sceneReFl);
    
    sensor = sensorCreate;
    sensor = sensorSetSizeToFOV(sensor,[sceneGet(sceneRe,'fov horizontal') sceneGet(sceneRe,'fov vertical')],sceneRe,oi);
    
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

    lRGBRe = ipGet(ipRe,'sensor channels');
    lRGBFl = ipGet(ipFl,'sensor channels');
    lRGBReFl = ipGet(ipReFl,'sensor channels');

    
    
    fName = fullfile(saveDir,sprintf('Light_%i_re.png',i));
    imwrite(lRGBRe,fName);
    fName = fullfile(saveDir,sprintf('Light_%i_fl.png',i));
    imwrite(lRGBFl,fName);
    fName = fullfile(saveDir,sprintf('Light_%i_reFl.png',i));
    imwrite(lRGBReFl,fName);
    
    figure; imshow([lRGBRe lRGBFl lRGBReFl]);
    
    vcAddObject(ipRe);
    vcAddObject(ipFl);
    vcAddObject(ipReFl);

    ipWindow;
    
end