close all;
clear variables;
clc;

inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

alpha = 0.1;
beta = 0.1;
nu = 0.1;

nNoiseLevels = 20;
nInstances = 10;
noiseLevels = logspace(-4,2,nNoiseLevels);



% Create basis function sets
nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = createBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = createBasisSet('excitation','wave',wave','n',nExBasis);
[emBasis, emScore] = createBasisSet('emission','wave',wave','n',nEmBasis);


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
   

nSamples = size(measVals,3);

reflEst = cell(nNoiseLevels,nInstances);
dMatEst = cell(nNoiseLevels,nInstances);
reflValsEst = cell(nNoiseLevels,nInstances);
flValsEst = cell(nNoiseLevels,nInstances);
measValsNoise = cell(noiseLevels,nInstances);

SNR = cell(nNoiseLevels,1);

try
    matlabpool open local
catch
end

for nl=1:nNoiseLevels
    
    SNR{nl} = measVals./noiseLevels(nl);
    
    for i=1:nInstances
        
        measValsNoise{nl,i} = max(measVals + randn(size(measVals))*noiseLevels(nl),0);
        
        localCameraGain = repmat(cameraGain,[1 1 nSamples]);
        localCameraOffset = repmat(cameraOffset,[1 1 nSamples]);
        
        nF = max(max(measValsNoise,[],1),[],2);
        localCameraGain = localCameraGain./repmat(nF,[nFilters nChannels 1]);
        measValsNoise = measValsNoise./repmat(nF,[nFilters nChannels 1]);
        
        [ reflEst{nl,i}, ~, ~, ~, ~, ~, dMatEst{nl,i}, reflValsEst{nl,i}, flValsEst{nl,i}, hist  ] = ...
            fiRecReflAndMultiFl( measValsNoise{nl,i}, camera, illuminant, localCameraGain*deltaL,...
            localCameraOffset, reflBasis, emBasis, exBasis, alpha, beta, beta, nu, 'maxIter',250);
        
    end
    
end

fName = fullfile(fiToolboxRootPath,'results','evaluation','multiFl_SNR.mat');
save(fName,'reflEst','dMatEst','reflValsEst','flValsEst','SNR',...
           'alpha','beta','nu','dMatRef','reflRef','inFName','measValsNoise');

try 
    matlabpool close
catch
end





