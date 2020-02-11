close all;
clear all;
clc;

ieInit;
[rootPath, parentPath] = fiToolboxRootPath();

wave = 400:5:700;


fName = fullfile(parentPath,'PhotosWithMacbeth','RAW','IMG_8431.CR2');

cmd = sprintf('/usr/local/bin/dcraw -r 1 1 1 1 -g 1 1 -v -d -4 -M -H 0 "%s"',fName);
system(cmd);

[path, name] = fileparts(fName);
[superPath, rawName] = fileparts(path);

if ~exist(fullfile(superPath,'Corrected'),'dir');
    mkdir(fullfile(superPath,'Corrected'));
end

switch name
    case 'IMG_8431'
        cp = [2089 3556;3390 3515;3363 2670;2068 2717];
        lightID = 4;
        % For narrowband light take the matrix from wide
        CM = [2.44045064547286 0.0384538073502225 0.0122113884354827;-0.00767895011651348 2.44216904577539 0.00169268181230433;0.0126554250438228 0.027009518827524 2.530954801508];
    case 'IMG_8422'
        cp = [2080 3554;3382 3540;3348 2688;2067 2716];
        lightID = 1;
    case 'IMG_8430'
        lightID = 3;
        cp = [2101 3567;3382 3554;3354 2709;2074 2736];
    case 'IMG_8427'
        lightID = 2;
        cp = [2087 3561;3389 3540;3341 2688;2067 2723];
end

rawFileName = fullfile(path,sprintf('%s.pgm',name));
rawImage = double(imread(rawFileName));
rawImage = rawImage/max(rawImage(:));

fName = fullfile(fiToolboxRootPath,'data','CanonG7x');
cfa = ieReadSpectra(fName,wave);

sensor = sensorCreate('bayer (gbrg)');
sensor = sensorSet(sensor,'wave',wave);
sensor = sensorSet(sensor,'filter transmissivities',cfa);
sensor = sensorSet(sensor,'size',[size(rawImage,1) size(rawImage,2)]);
sensor = sensorSet(sensor,'volts',rawImage);
sensor = sensorSet(sensor,'noise flag',0);

ieAddObject(sensor);
sensorWindow();

[realMacbethPatch, ~, ~, cp] = macbethSelect(sensor,1,1, cp);
realMacbethPatch = cellfun(@nanmean,realMacbethPatch,'UniformOutput',false)';
realMacbethPatch = cell2mat(realMacbethPatch);

disp(mat2str(cp));

realIp = ipCreate();
realIp = ipSet(realIp,'correction method illuminant','none');
realIp = ipSet(realIp,'conversion method sensor','none');
realIp = ipCompute(realIp,sensor);

demosRealRawImage = ipGet(realIp,'sensor channels');
sz = size(demosRealRawImage);
demosRealRawImageVec = reshape(demosRealRawImage,[sz(1)*sz(2), sz(3)]); 

ieAddObject(realIp);
ipWindow();

%% Simulated data

path = fullfile(parentPath,'Photos');
fName = fullfile(path,sprintf('%s_RenderedLight_%i_Macbeth.png','Target6',lightID));
simMacbeth = double(imread(fName));
simMacbeth = simMacbeth/255;

ip = ipCreate;
ip.data.result = double(simMacbeth);
cp = [2, 94; 128, 94; 128, 11; 1, 11];

simMacbethPatch = macbethSelect(ip,1,1,cp);
simMacbethPatch = cellfun(@nanmean,simMacbethPatch,'UniformOutput',false)';
simMacbethPatch = cell2mat(simMacbethPatch);

if ~exist('CM','var')
    CM = realMacbethPatch\simMacbethPatch;
end
corrRealMacbethPatch = realMacbethPatch*CM;

figure;
hold on; grid on; box on;
plot(simMacbethPatch,corrRealMacbethPatch,'o');
xlabel('Simulated');
ylabel('Corrected measured');

realCorrVec = min(max(demosRealRawImageVec*CM,0),1);
realCorr = reshape(realCorrVec,sz);
figure; imshow(realCorr);

imwrite(realCorr,fullfile(superPath,'Corrected',sprintf('%s.png',name)));


