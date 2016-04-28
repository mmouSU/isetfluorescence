function [ reflectanceEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, donaldsonMat, predRefl, predFl, hist  ] = fiRecReflAndMultiFl( measVals, cameraMat, illuminant, cameraGain, cameraOffset, basisRefl, basisEm, basisEx, alpha, beta, gamma, eta, varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('dMatRef',[]);
p.addParamValue('reflRef',[]);
p.parse(varargin{:});
inputs = p.Results;


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
    if ~isempty(inputs.dMatRef)
        dMatRef = inputs.dMatRef{i};
    else
        dMatRef = [];
    end
    if ~isempty(inputs.reflRef);
        reflRef = inputs.reflRef(:,i);
    else
        reflRef = [];
    end

    [reflectanceEst(:,i), rfCoeffs(:,i), emEst(:,i), emCoeffs(:,i), exEst(:,i), exCoeffs(:,i), donaldsonMat{i}, predRefl(:,:,i), predFl(:,:,i), hist{i}] = ...
    fiRecOneReflAndMultiFl(input,cameraMat,illuminant,cameraGain(:,:,i),basisRefl,basisEm,basisEx,alpha,beta,gamma,eta,varargin{:},'dMatRef',dMatRef,'reflRef',reflRef);
    
    fprintf('Done!\n');
end

end

