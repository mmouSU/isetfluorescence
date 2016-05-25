% Display images for all pairs of camera filter and illuminant captured for 
% a particular scene. This script generates images in Fig. 11 in the paper.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
% saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');
saveDir = [];


testFileName = 'Macbeth+Fl';
backgroundFileName = 'Background';

wave = 380:4:1000;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Load the light spectra (in photons). We scale by the maximum for
% numerical stability. This does not matter in the long run as we are
% calibrating for an unknown gain paramter.
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
illuminantPhotons = Energy2Quanta(wave,illuminant);
illuminantPhotons = illuminantPhotons/max(illuminantPhotons(:));

nChannels = size(illuminant,2);

% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);



%% Calibration
% The whole linear image formation model has an unknown gain parameter, we
% are computing this gain for every pixel by taking an image of a white
% surface with a known reflectance. The gain is just the result of the
% division between pixel intensities and linear model predictions.

% Predict the response of a target surface
prediction = deltaL*((camera')*illuminantPhotons);

% Generate the gain map for every pixel
fName = fullfile(fiToolboxRootPath,'data','experiments',backgroundFileName);
[~, ~, scaledRAW, shutterBackground] = fiReadImageStack(fName);
hh = size(scaledRAW,1);
ww = size(scaledRAW,2);

% For every illuminant channel and the monochromatic filter pick a
% reference point (say in the middle). Compute how would you need to scale
% other pixels so that the image intensity is uniform. We are using the
% monochromatic channel for the computations and then using the same map
% across all the channels
refPt = scaledRAW(hh/2,ww/2,1,:);
refImg = repmat(refPt,[hh ww 1 1]);
scaleMap = refImg./scaledRAW(:,:,1,:);
scaleMap = repmat(scaleMap,[1 1 nFilters 1]);

nonuniformityCorrected = scaledRAW.*scaleMap;


%% Extract data from a Macbeth image
fName = fullfile(fiToolboxRootPath,'data','experiments',testFileName);
[~, ~, scaledMacbeth, shutterMacbeth] = fiReadImageStack(fName);
linearVals = scaledMacbeth.*scaleMap;

for f=1:nFilters
    
    nF = linearVals(:,:,f,:);
    nF = max(nF(:));
    
    for c=1:nChannels
        fi = figure;
        tmp = max(imresize(linearVals(:,:,f,c),0.2),0);
        imshow((tmp/nF).^(1/2.2),'InitialMagnification',100,'Border','tight');
        
        if ~isempty(saveDir)
            fName = fullfile(saveDir,sprintf('image_f%i_c%i.eps',f,c));
            print('-depsc',fName);
            close(fi);
        end
        
    end
end


