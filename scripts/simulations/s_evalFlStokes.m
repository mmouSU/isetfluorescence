close all;
clear all;
clc;

close all;
clear all;
clc;

% Evaluate the accuracy on the entire dataset

dataset = 'McNamara-Boswell';
wave = 380:4:1000;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

alpha = 0.1;
beta = 0.1;
height = 4;
width = 6;
nSamples = height*width;
nFluorophores = 1;

flQe = 0.1;
stokesRanges = [0 25 50 75 100];

nStokes = length(stokesRanges)-1;

% Create basis function sets
nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = createBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = createBasisSet('excitation','wave',wave','n',nExBasis);
[emBasis, emScore] = createBasisSet('emission','wave',wave','n',nEmBasis);


totalPixelErr = zeros(nStokes,1);
reflPixelErr = zeros(nStokes,1);
flPixelErr = zeros(nStokes,1);
reflErr = zeros(nStokes,1);
exErr = zeros(nStokes,1);
exNormErr = zeros(nStokes,1);
emErr = zeros(nStokes,1);
emNormErr = zeros(nStokes,1);

totalPixelStd = zeros(nStokes,1);
reflPixelStd = zeros(nStokes,1);
flPixelStd = zeros(nStokes,1);
reflStd = zeros(nStokes,1);
exStd = zeros(nStokes,1);
exNormStd = zeros(nStokes,1);
emStd = zeros(nStokes,1);
emNormStd = zeros(nStokes,1);

%% The main cross-validation loop

try
    matlabpool open local
catch 
end

parfor i=1:nStokes

    fName = fullfile(fiToolboxRootPath,'data','simulations',sprintf('%s_%ix%ix%i_stokes_%i.mat',dataset,height,width,nFluorophores,stokesRanges(i+1)));
    data = load(fName);

    
    cameraGain = repmat(data.cameraGain,[1 1 nSamples]);
    cameraOffset = repmat(data.cameraOffset,[1 1 nSamples]);
    
    
    [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndFl( data.measVals, data.camera, cameraGain*deltaL, cameraOffset, data.illuminantPhotons,...
        reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter', 250 );

    
    %% Evaluation
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [totalPixelErr(i), totalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(data.measVals,[nChannels*nFilters,nSamples]), '');
    [reflPixelErr(i), reflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(data.reflValsRef,[nChannels*nFilters,nSamples]), '');
    [flPixelErr(i), flPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(data.flValsRef,[nChannels*nFilters,nSamples]), '');
    
    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, data.reflRef, '');
    
    [emErr(i), emStd(i)] = fiComputeError(emEst, data.emRef, '');
    [emNormErr(i), emNormStd(i)] = fiComputeError(emEst, data.emRef, 'normalized');
    
    [exErr(i), exStd(i)] = fiComputeError(exEst, data.exRef, '');
    [exNormErr(i), exNormStd(i)] = fiComputeError(exEst, data.exRef, 'normalized');
    

end

try
    matlabpool close
catch
end

%% Save results
dirName = fullfile(fiToolboxRootPath,'results','evaluation');
if ~exist(dirName,'dir'), mkdir(dirName); end

fName = fullfile(dirName,[dataset '_simStokes_Fl.mat']);

save(fName,'stokesRanges',...
           'totalPixelErr','reflPixelErr','flPixelErr','reflErr','exErr','exNormErr','emErr','emNormErr',...
           'totalPixelStd','reflPixelStd','flPixelStd','reflStd','exStd','exNormStd','emStd','emNormStd');

       
