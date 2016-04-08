close all;
clear all;
clc;

inFName = 'McNamara-Boswell_4x6x1_qe_0.10';

% Load simulation data
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

alpha = 0.01;
sigma = 0.001;

% Create basis function sets
nReflBasis = 5;
nExBasis = 12;

[reflBasis, reflScore] = createBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = createBasisSet('excitation','wave',wave','n',nExBasis);


% Load the light spectra (in photons)
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
illuminant = Energy2Quanta(wave,illuminant);
nChannels = size(illuminant,2);

% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);


% Re-scale the measured values so that the maximum for each measurement is
% 1
       
nSamples = size(measVals,3);
cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

nF = max(max(measVals,[],1),[],2);
measVals = measVals./repmat(nF,[nFilters nChannels 1]);
cameraGain = cameraGain./repmat(nF,[nFilters nChannels 1]);
%%

[ reflEst, emEst, exEst, dMatEst, reflValsEst, flValsEst, hist ] = fiRecReflAndFlNucNorm( measVals,...
    camera, cameraGain*deltaL, cameraOffset, illuminant, alpha, sigma, 'maxIter',100 );

measValsEst = reflValsEst + flValsEst;


[err, std] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals - cameraOffset,[nChannels*nFilters,nSamples]), '');
fprintf('Total pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), '');
fprintf('Reflected pixel error %.3f, std %.3f\n',err,std);

[err, std] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), '');
fprintf('Fluoresced pixel error %.3f, std %.3f\n',err,std);

% [err, std] = fiComputeError(reflEst, reflRef, '');
% fprintf('Reflectance error %.3f, std %.3f\n',err,std);





%% Plot the results

fName = fullfile(fiToolboxRootPath,'data','flCmap');
load(fName);

% Predicted vs. simulated pixel intensities

figure;
hold all; grid on; box on;
plot(measValsEst(:),measVals(:),'.');
xlabel('Model predicted pixel value');
ylabel('ISET pixel value');

figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on; axis square;

    tmp1 = measValsEst(:,:,sampleID);
    tmp2 = measVals(:,:,sampleID);
    plot(tmp1(:),tmp2(:),'.');

end
end


% Convergence
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot([hist{sampleID}.prRes, hist{sampleID}.dualRes],'LineWidth',2);
    

end
end


% Estimated vs. ground truth reflectance
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot(wave,diag(reflEst{sampleID}),'g','LineWidth',2);
    plot(wave,reflRef(:,sampleID),'b--','LineWidth',2);
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);

    rmse = sqrt(mean((diag(reflEst{sampleID}) - reflRef(:,sampleID)).^2));
    title(sprintf('RMSE %.2f',rmse));

end
end

% Estimated vs. ground truth Donaldson matrices: scatter plot
figure;
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;

    plot(dMatEst{sampleID}(:),dMatRef{sampleID}(:),'.');

    rmse = sqrt(mean((dMatEst{sampleID}(:) - dMatRef{sampleID}(:)).^2));
    title(sprintf('RMSE %.2e',rmse));
end
end

% Estimated vs. ground truth Donaldson matrics: scale
figure;
set(gcf,'Colormap',flCmap);
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    
    data = [dMatEst{sampleID} dMatRef{sampleID}];
    imagesc(data);

end
end


% Estimated vs. ground truth Donaldson matrics: shape
figure;
set(gcf,'Colormap',flCmap);
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    
    data = [dMatEst{sampleID}/max(dMatEst{sampleID}(:)) dMatRef{sampleID}/max(dMatRef{sampleID}(:))];
    imagesc(data);

end
end


