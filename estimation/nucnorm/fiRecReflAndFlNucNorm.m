function [ reflEst, emEst, exEst, dMatEst, predRefl, predFl, hist  ] = fiRecReflAndFlNucNorm( measVals, cameraMat, cameraGain, cameraOffset, illuminant, alpha, sigma, varargin )

nSamples = size(measVals,3);
nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);
nWaves = size(cameraMat,1);

reflEst = cell(nSamples,1);

emEst = zeros(nWaves,nSamples);
exEst = zeros(nWaves,nSamples);
dMatEst = cell(nSamples,1);

predRefl = zeros([nFilters, nChannels, nSamples]);
predFl = zeros([nFilters, nChannels, nSamples]);

hist = cell(nSamples,1);

for i=1:nSamples

    input = measVals(:,:,i) - cameraOffset(:,:,i);
    
    [reflEst{i}, emEst(:,i), exEst(:,i), dMatEst{i}, predRefl(:,:,i), predFl(:,:,i), hist{i}] = ...
    fiRecOneReflAndFlNucNorm(input,cameraMat,cameraGain(:,:,i),illuminant,alpha,sigma,varargin{:});
    
end

end

