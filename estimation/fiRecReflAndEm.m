function [ reflEst, rfCoeffs, emEst, emCoeffs, emWghts, predRefl, predFl, hist  ] = fiRecReflAndEm( measVals, cameraMat, cameraGain, cameraOffset, illuminant, basisRefl, basisEm, alpha, beta, varargin )

nSamples = size(measVals,3);
nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);
nWaves = size(cameraMat,1);
nBasisRefl = size(basisRefl,2);
nBasisEm = size(basisEm,2);

reflEst = zeros(nWaves,nSamples);
rfCoeffs = zeros(nBasisRefl,nSamples);

emEst = zeros(nWaves,nSamples);
emCoeffs = zeros(nBasisEm,nSamples);
emWghts = zeros(nChannels,nSamples);

predRefl = zeros([nFilters, nChannels, nSamples]);
predFl = zeros([nFilters, nChannels, nSamples]);

hist = cell(nSamples,1);

for i=1:nSamples
    
    fprintf('Processing sample %i/%i ... ',i,nSamples);
    inputs = measVals(:,:,i) - cameraOffset(:,:,i);

    [reflEst(:,i), rfCoeffs(:,i), emEst(:,i), emCoeffs(:,i), emWghts(:,i), predRefl(:,:,i), predFl(:,:,i), hist{i}] = ...
    fiRecOneReflAndEm(inputs,cameraMat,cameraGain(:,:,i),illuminant,basisRefl,basisEm,alpha,beta,varargin{:});

    nIter = length(hist{i}.objValsReEm);

    fprintf('Done (%i iterations)\n',nIter);

end

end

