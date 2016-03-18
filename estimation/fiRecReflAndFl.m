function [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, predRefl, predFl, hist  ] = fiRecReflAndFl( measVals, cameraMat, cameraGain, cameraOffset, illuminant, basisRefl, basisEm, basisEx, alpha, beta, gamma, varargin )

nSamples = size(measVals,3);
nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);
nWaves = size(cameraMat,1);
nBasisRefl = size(basisRefl,2);
nBasisEm = size(basisEm,2);
nBasisEx = size(basisEx,2);

reflEst = zeros(nWaves,nSamples);
rfCoeffs = zeros(nBasisRefl,nSamples);

emEst = zeros(nWaves,nSamples);
emCoeffs = zeros(nBasisEm,nSamples);

exEst = zeros(nWaves,nSamples);
exCoeffs = zeros(nBasisEx,nSamples);

predRefl = zeros([nFilters, nChannels, nSamples]);
predFl = zeros([nFilters, nChannels, nSamples]);

hist = cell(nSamples,1);

for i=1:nSamples
    fprintf('Processing sample %i/%i ... ',i,nSamples);

    input = measVals(:,:,i) - cameraOffset(:,:,i);
    
    [reflEst(:,i), rfCoeffs(:,i), emEst(:,i), emCoeffs(:,i), exEst(:,i), exCoeffs(:,i), predRefl(:,:,i), predFl(:,:,i), hist{i}] = ...
        fiRecOneReflAndFl(input,cameraMat,cameraGain(:,:,i),illuminant,basisRefl,basisEm,basisEx,alpha,beta,gamma,varargin{:});

    nIter = length(hist{i}.objValsReEm);

    fprintf('Done (%i iterations)\n',nIter);

end


end

