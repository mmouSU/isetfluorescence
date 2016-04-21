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

flQe = logspace(-3,0,10);
nQe = length(flQe);

% Create basis function sets
nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = createBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = createBasisSet('excitation','wave',wave','n',nExBasis);
[emBasis, emScore] = createBasisSet('emission','wave',wave','n',nEmBasis);


totalPixelErr = zeros(nQe,1);
reflPixelErr = zeros(nQe,1);
flPixelErr = zeros(nQe,1);
reflErr = zeros(nQe,1);
exErr = zeros(nQe,1);
exNormErr = zeros(nQe,1);
emErr = zeros(nQe,1);
emNormErr = zeros(nQe,1);

totalPixelStd = zeros(nQe,1);
reflPixelStd = zeros(nQe,1);
flPixelStd = zeros(nQe,1);
reflStd = zeros(nQe,1);
exStd = zeros(nQe,1);
exNormStd = zeros(nQe,1);
emStd = zeros(nQe,1);
emNormStd = zeros(nQe,1);

%% The main cross-validation loop

try
    matlabpool open local
catch 
end

parfor i=1:nQe

    fName = fullfile(fiToolboxRootPath,'data','simulations',sprintf('%s_%ix%ix%i_qe_%.3f.mat',dataset,height,width,nFluorophores,flQe(i)));
    data = load(fName);

    
    cameraGain = repmat(data.cameraGain,[1 1 nSamples]);
    cameraOffset = repmat(data.cameraOffset,[1 1 nSamples]);
    
    
    [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndFl( data.measVals, data.camera, cameraGain*deltaL, cameraOffset, data.illuminantPhotons,...
        reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter',25 );

    
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

fName = fullfile(dirName,[dataset '_simQe_Fl.mat']);

save(fName,'flQe',...
           'totalPixelErr','reflPixelErr','flPixelErr','reflErr','exErr','exNormErr','emErr','emNormErr',...
           'totalPixelStd','reflPixelStd','flPixelStd','reflStd','exStd','exNormStd','emStd','emNormStd');

       
