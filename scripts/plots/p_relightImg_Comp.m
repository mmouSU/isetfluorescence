% Compare images of a test targed under different illuminants. This script
% generates synthetic images in Fig. 6 in the Supplemental Material.
%
% Note, you should've generated synthetic data running
% p_relightImg_Natural.m
%
% Copyright, Henryk Blasinski

close all;
clear all;
clc;

ieInit;
wave = 380:4:1000; wave = wave(:);
fs = 6;

targetName = 'Natural';

% Define path to directory with the outputs of the p_relightImg_Natural.m
simPath = fullfile('~','Desktop','Test images');

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
saveDir = fullfile('~','Desktop','Natural-Corrected');
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

switch targetName
    case 'Natural'

        % Define path to directory with photos captured with Canon G7x
        capPath = fullfile(fiToolboxRootPath,'..','Photos','Natural');
        
        % Load illuminants
        % These are just placeholder illuminants.
        fName = fullfile(fiToolboxRootPath,'camera','illuminants');
        illuminant = ieReadSpectra(fName,wave);
        illSubset = illuminant(:,[2 6]);

        % Images under these illuminants are included in the Supplemental Materials
        fName = fullfile(fiToolboxRootPath,'data','Broadband9800K');
        illSubset = [illSubset, ieReadSpectra(fName,wave)]; 
        fName = fullfile(fiToolboxRootPath,'data','Broadband6500K');
        illSubset = [illSubset, ieReadSpectra(fName,wave)];
        fName = fullfile(fiToolboxRootPath,'data','Broadband2000K');
        illSubset = [illSubset, ieReadSpectra(fName,wave)];
        fName = fullfile(fiToolboxRootPath,'data','NarrowbandBlue');
        illSubset = [illSubset, ieReadSpectra(fName,wave)];

        nIlls = size(illSubset,2);



        % Corner points for Macbth chart images
        % Captured:
        cp = cell(nIlls,1);
        cp{1} = [2785 2081;3400 1506;2982 1063;2383 1654];
        cp{2} = [2760 2048;3384 1490;2990 1071;2366 1638];
        cp{3} = [2769 2015;3359 1457;2998 1055;2375 1630];
        cp{4} = [2727 2040;3343 1474;2965 1088;2350 1654];
        cp{5} = [2777 2064;3384 1498;2990 1096;2383 1679];
        cp{6} = [2670 1982;3294 1441;2933 1039;2317 1588];

        % Simulated:
        simCp = [2 178;241 178;240 18;2 19];

end

% Load Macbeth chart reflectance
fName = fullfile(fiToolboxRootPath,'data','macbethChart');
refl = ieReadSpectra(fName,wave);

err = cell(nIlls,1);
errRe = cell(nIlls,1);

measMacbeth = zeros(3,24*nIlls);
simMacbeth = zeros(3,24*nIlls);

for i=3:nIlls
    
    %% Get the Canon image of the Macbeth chart
    macbeth = double(imread(fullfile(capPath,sprintf('CapturedLight_%i_Macbeth.tiff',i))))/255;
    
    sensor = sensorCreate('monochrome');
    sensor = sensorSet(sensor,'size',[size(macbeth,1) size(macbeth,2)]);
    
    for j=1:3
        sensor = sensorSet(sensor,'volts',macbeth(:,:,j));
        ieAddObject(sensor);
        sensorWindow;
        [vals,~,~,cps] = macbethSelect(sensor,1,1,cp{i});
        measMacbeth(j,(i-1)*24+1:i*24) = cellfun(@nanmean,vals);
    end
    
    %% Get the simulated image of the Macbeth chart
    macbeth = double(imread(fullfile(simPath,sprintf('%s_renderedLight_%i_Macbeth.png',lower(targetName),i))))/255;
    
    sensor = sensorCreate('monochrome');
    sensor = sensorSet(sensor,'size',[size(macbeth,1) size(macbeth,2)]);
    
    for j=1:3
        sensor = sensorSet(sensor,'volts',macbeth(:,:,j));
        ieAddObject(sensor);
        sensorWindow;
        [vals,~,~,cps] = macbethSelect(sensor,0,1,simCp);
        simMacbeth(j,(i-1)*24+1:i*24) = cellfun(@nanmean,vals);
    end
end
    
% Least-squares fit of a 3x3 color conversion matrix.
sim2canon = simMacbeth'\measMacbeth';

tmp = (simMacbeth'*sim2canon)';

figure; 
hold on; grid on; box on;
plot(measMacbeth(:),tmp(:),'.');
xlabel('Measured pixel values');
ylabel('Corrected pixel values');
    

%%
% For captures under each illuminant apply a 3x3 conversion computed from
% the images of a macbeth chart and compute pixelwise error maps. 
for i = 3:nIlls
    
    cImg = im2double(imread(fullfile(capPath, sprintf('CapturedLight_%i_cropped.png',i))));
    xImg = im2double(imread(fullfile(simPath, sprintf('%s_renderedLight_%i_reFl.png',lower(targetName),i))));
    xImgFl = im2double(imread(fullfile(simPath, sprintf('%s_renderedLight_%i_fl.png',lower(targetName),i))));
    xImgRe = im2double(imread(fullfile(simPath, sprintf('%s_renderedLight_%i_re.png',lower(targetName),i))));

    %Simulate bilinear demosaicing
    R = cImg(:,:,1);
    R(1:2:end,1:2:end) = 0;
    G1 = cImg(:,:,2);
    G1(2:2:end,1:2:end) = 0;
    G2 = cImg(:,:,2);
    G2(1:2:end,2:2:end) = 0;
    B = cImg(:,:,3);
    B(2:2:end,2:2:end) = 0;
 
    RGB = cat(3,R,G1+G2,B);
 
    demosaiced = ieBilinear(RGB,[1 2; 2 3]);
    demosaiced = demosaiced/max(demosaiced(:));
    cImg = imfilter(demosaiced,fspecial('gaussian',[11 11],1.25));
    cImg = min(max(cImg,0),1);
    
    % Transform
    xImgRe = imageLinearTransform(xImgRe,sim2canon);
    xImgFl = imageLinearTransform(xImgFl,sim2canon);
    xImg = imageLinearTransform(xImg,sim2canon);

    % Gain
    scale = xImg(:)\cImg(:);
    xImg = max(xImg*scale, 0);
    xImgRe = max(xImgRe*scale, 0);
    xImgFl = max(xImgFl*scale, 0);

    figure; imshow(xImgRe.^(1/2.2),'Border','tight');
    if ~isempty(saveDir) 
        print('-depsc',fullfile(saveDir, sprintf('%s_renderedLight_%i_re.eps',lower(targetName),i))); 
    end
    
    figure; imshow(xImgFl.^(1/2.2),'Border','tight');
    if ~isempty(saveDir)
        print('-depsc',fullfile(saveDir, sprintf('%s_renderedLight_%i_fl.eps',lower(targetName),i))); 
    end
    
    figure; imshow(xImg.^(1/2.2),'Border','tight');
    if ~isempty(saveDir)
        print('-depsc',fullfile(saveDir, sprintf('%s_renderedLight_%i_reFl.eps',lower(targetName),i))); 
    end

    figure; imshow(cImg.^(1/2.2),'Border','tight');
    if ~isempty(saveDir)
        print('-depsc',fullfile(saveDir, sprintf('%s_capturedLight_%i.eps',lower(targetName),i))); 
    end
    
    figure; imshow([xImgRe xImgFl xImg cImg].^(1/2.2),'Border','tight');

    err{i} = sqrt(nanmean((cImg - xImg).^2,3));
    errRe{i} = sqrt(nanmean((cImg - xImgRe).^2,3));
    figure; imagesc(err{i}); colorbar;
end

%% Display error maps

h = size(err{i},1);
w = size(err{i},2);

minVal = 0;
maxVal = 1;
load('flCmap.mat');

for i=3:6
    
    figure;
    imshow(err{i},[minVal maxVal],'Border','tight');
    colormap(flCmap);
    
    fprintf('Average camera rgb %f, median camera rgb %f, mode camera rgb %f\n',mean(err{i}(:)),median(err{i}(:)),mode(err{i}(:)));
    
    
    axis image;
    axis off;
    set(gca,'fontsize',6);
    set(gcf,'PaperUnits','Centimeters');
    set(gcf,'PaperPosition',[1 1 6*[1 h/w]]);
    
    if ~isempty(saveDir)
        print('-depsc',fullfile(saveDir, sprintf('%s_light_%i_error.eps',lower(targetName),i))); 
    end
    
    sz = get(gcf,'PaperPosition');
    
    
    % Colorbar
    if i==6
        figure;
        set(gcf,'Colormap',flCmap);
        
        left=100; bottom=200 ; width=50 ; height=500;
        pos=[left bottom width height];
        axis off
        
        
        cb = colorbar('position',[0.1 0.15  0.2  0.8]);
        set(gca,'CLim',[minVal maxVal]);
        set(gca,'fontsize',1.5*fs);
        set(cb,'TickLabelInterpreter','latex');
        set(cb,'colormap',flCmap);
        set(gcf,'OuterPosition',pos);
        set(gcf,'Units','centimeters');
        set(gcf,'PaperPosition',[1 1 2 sz(4)]);
        if ~isempty(saveDir)
            print('-depsc',fullfile(saveDir, sprintf('errorbar.eps')));
        end
    end
end



