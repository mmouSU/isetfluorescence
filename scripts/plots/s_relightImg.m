close all;
clear all;
clc;

ieInit;
wave = 380:4:1000;
nWaves = length(wave);

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures','Canon G7x V2');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

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


fName = fullfile(fiToolboxRootPath,'data','macbethChart');
refl = ieReadSpectra(fName,wave);
sceneMacbeth = sceneReflectanceArray(reshape(refl',[4 6 nWaves]),10,wave);

% Camera filter transmissivities
fName = fullfile(isetRootPath,'data','sensor','colorfilters','Canon600D');
cfa = ieReadSpectra(fName,wave);


%%

for i=1:size(illSubset,2)      
        
    sceneMacbeth = sceneAdjustIlluminant(sceneMacbeth,illSubset(:,i),0);
    vcAddObject(sceneMacbeth);
    
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
    
    
    fName = fullfile(saveDir,sprintf('RenderedLight_%i_re.png',i));
    imwrite(lRGBRe,fName);
    fName = fullfile(saveDir,sprintf('RenderedLight_%i_fl.png',i));
    imwrite(lRGBFl,fName);
    fName = fullfile(saveDir,sprintf('RenderedLight_%i_reFl.png',i));
    imwrite(lRGBReFl,fName);
    fName = fullfile(saveDir,sprintf('RenderedLight_%i_Macbeth.png',i));
    imwrite(lRGBMacbeth,fName);
    
    fName = fullfile(saveDir,sprintf('RenderedLight_%i_re_XYZ.mat',i));
    XYZre = imresize(ieXYZFromPhotons(sceneGet(sceneRe,'photons'),wave),[size(lRGBRe,1) size(lRGBRe,2)]);
    save(fName,'XYZre');
    fName = fullfile(saveDir,sprintf('RenderedLight_%i_fl_XYZ.mat',i));
    XYZfl = imresize(ieXYZFromPhotons(sceneGet(sceneFl,'photons'),wave),[size(lRGBFl,1) size(lRGBFl,2)]);
    save(fName,'XYZfl');
    fName = fullfile(saveDir,sprintf('RenderedLight_%i_reFl_XYZ.mat',i));
    XYZreFl = imresize(ieXYZFromPhotons(sceneGet(sceneReFl,'photons'),wave),[size(lRGBReFl,1) size(lRGBReFl,2)]);
    save(fName,'XYZreFl');

    
    
    figure; imshow([lRGBRe lRGBFl lRGBReFl lRGBMacbeth]);
    figure; imshow([sRGBRe sRGBFl sRGBReFl sRGBMacbeth]);

    
    vcAddObject(ipRe);
    vcAddObject(ipFl);
    vcAddObject(ipReFl);

    ipWindow;
    drawnow;
end