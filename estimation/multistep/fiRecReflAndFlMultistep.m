function [ reflEst, reflCoeffs, emEst, emChromaticity, exEst, exCoeffs, predRefl, predFl, hist ] = fiRecReflAndFlMultistep( measVals, cameraMat, cameraGain, cameraOffset, illuminant, basisRefl, DB, basisEx, alpha, gamma )

nSamples = size(measVals,3);
nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);
nWaves = size(cameraMat,1);
nBasisRefl = size(basisRefl,2);
nBasisEx = size(basisEx,2);

reflEst = zeros(nWaves,nSamples);
reflCoeffs = zeros(nBasisRefl,nSamples);

emChromaticity = zeros(nFilters,nSamples);

emEst = zeros(nWaves,nSamples);

exEst = zeros(nWaves,nSamples);
exCoeffs = zeros(nBasisEx,nSamples);

predRefl = zeros([nFilters, nChannels, nSamples]);
predFl = zeros([nFilters, nChannels, nSamples]);

hist = cell(nSamples,1);

% Scale to avoid numerical issues
scaleFac = max(illuminant(:));
illuminant = illuminant/scaleFac;
cameraGain = cameraGain*scaleFac;


for i=1:nSamples
    
    fprintf('Processing sample %i ... ',i);
    
    input = measVals(:,:,i) - cameraOffset(:,:,i);
    
    [reflEst(:,i), reflCoeffs(:,i), emChromaticity(:,i), hist{i}] = fiRecOneRefl(input,...
        cameraMat,cameraGain(:,:,i),illuminant,basisRefl,alpha);
    
    
    [exEst(:,i), exCoeffs(:,i)] = fiRecOneEx(input,reflCoeffs(:,i),emChromaticity(:,i),...
        cameraMat,cameraGain(:,:,i),illuminant,basisRefl,basisEx,gamma);
    
    [emEst(:,i)] = fiRecOneEm(emChromaticity(:,i),cameraMat,DB);
    
    predRefl(:,:,i) = cameraGain(:,:,i).*(cameraMat'*diag(reflEst(:,i))*illuminant);
    
    flEx = exEst(:,i)'*illuminant;
    flEm = cameraMat'*emEst(:,i);
    
    % We are trying to find a single gain parameter on fluorescence emission
    % that best explains the data. This parameter is roughly equivalent to
    % practical quantum efficiency.
    cvx_begin quiet
        variable wf
        minimize norm(input - predRefl(:,:,i) - cameraGain(:,:,i).*(flEm*flEx)*wf,'fro')
        subject to
            wf >= 0
    cvx_end

    fprintf(' (%i iterations) Done!\n',length(hist{i}.objValsRefl));
    
    predFl(:,:,i) = cameraGain(:,:,i).*(flEm*flEx).*wf;
    
    nF = max(exEst(:,i));
    exEst(:,i) = exEst(:,i)/nF;
    emEst(:,i) = emEst(:,i)*nF*wf;
    
    % tmp = cameraGain(:,:,i).*(cameraMat'*(emEst(:,i)*exEst(:,i)')*illuminant)
end


end