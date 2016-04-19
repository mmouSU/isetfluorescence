close all;
clear all;
clc;

ieInit;

wave = 380:4:1000;

fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);

illSubset = illuminant(:,[1 3]);
tmp = illuminantCreate('blackbody',wave,10000);
illSubset = [illSubset illuminantGet(tmp,'energy')];
tmp = illuminantCreate('blackbody',wave,6500);
illSubset = [illSubset illuminantGet(tmp,'energy')];
tmp = illuminantCreate('blackbody',wave,2000);
illSubset = [illSubset illuminantGet(tmp,'energy')];


dataDir = fullfile(fiToolboxRootPath,'results','experiments');
c = 128;
r = 103; 

% Load scene reflectance
reflArray = zeros(r,c,156);
for cc=1:c
    fName = fullfile(dataDir,sprintf('multiFl_Natural_col_%i.mat',cc));
    data = load(fName);
    reflArray(:,cc,:) = data.reflEst';
    
    for rr=1:r
        deltaL = data.wave(2) - data.wave(1);
        flArray(rr,cc) = fluorophoreCreate('type','fromdonaldsonmatrix',...
                                         'DonaldsonMatrix',data.dMatEst{rr}/deltaL,...
                                         'wave',data.wave);
    end
    
end    

%% Create a scene with the specified reflectance properties and a particular
% illuminant
sceneReTemplate = sceneReflectanceArray(reflArray,1,wave);
sceneFlTemplate = fluorescentSceneCreate('type','fromfluorophore','fluorophore',flArray);

%%

for i=1:1 %size(illSubset,2)      
        
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
    
    figure; imshow(sceneGet(sceneFl,'RGB'),'Border','tight');
    
    oi = oiCreate;
    oiRe = oiCompute(oi,sceneRe);
    oiFl = oiCompute(oi,sceneFl);
    oiReFl = oiCompute(oi,sceneReFl);
    
    sensor = sensorCreate;
    sensor = sensorSetSizeToFOV(sensor,9,sceneRe,oiRe);
    
    sensorRe = sensorCompute(sensor,oiRe);
    sensorFl = sensorCompute(sensor,oiFl);
    sensorReFl = sensorCompute(sensor,oiReFl);
    
    ip = ipCreate;
    ipRe = ipCompute(ip,sensorRe);
    ipRe = ipSet(ipRe,'name',sprintf('refl - ill%i',i));
    
    ipFl = ipCompute(ip,sensorFl);
    ipFl = ipSet(ipFl,'name',sprintf('fl - ill%i',i));
    ipReFl = ipCompute(ip,sensorReFl);
    ipReFl = ipSet(ipReFl,'name',sprintf('refl+fl - ill%i',i));

    
    figure; imshow(ipGet(ipReFl,'sensor channels').^(1/2.2));
    
    vcAddObject(ipRe);
    vcAddObject(ipFl);
    vcAddObject(ipReFl);

    ipWindow;
    
end