function [ reflectanceEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, donaldsonMat, predRefl, predFl, hist  ] = fiRecReflAndMultiFl( measVals, cameraMat, illuminant, cameraGain, cameraOffset, basisRefl, basisEm, basisEx, alpha, beta, gamma, eta, varargin)

nSamples = size(measVals,3);
nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);
nWaves = size(cameraMat,1);
nBasisRefl = size(basisRefl,2);
nBasisEm = size(basisEm,2);
nBasisEx = size(basisEx,2);

reflectanceEst = zeros(nWaves,nSamples);
rfCoeffs = zeros(nBasisRefl,nSamples);

emEst = zeros(nWaves,nSamples);
emCoeffs = zeros(nBasisEm,nSamples);

exEst = zeros(nWaves,nSamples);
exCoeffs = zeros(nBasisEx,nSamples);

donaldsonMat = cell(nSamples,1);

predRefl = zeros([nFilters, nChannels, nSamples]);
predFl = zeros([nFilters, nChannels, nSamples]);

hist = cell(nSamples,1);

for i=1:nSamples
    
    fprintf('Processing sample %i of %i\n',i,nSamples);
    
    input = measVals(:,:,i) - cameraOffset(:,:,i);

    [reflectanceEst(:,i), rfCoeffs(:,i), emEst(:,i), emCoeffs(:,i), exEst(:,i), exCoeffs(:,i), donaldsonMat{i}, predRefl(:,:,i), predFl(:,:,i), hist{i}] = ...
    fiRecOneReflAndMultiFl(input,cameraMat,illuminant,cameraGain(:,:,i),basisRefl,basisEm,basisEx,alpha,beta,gamma,eta,varargin{:});
    
    fprintf('Done!\n');
end

end

