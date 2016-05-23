function [ reflEst, reflCoeffs, emEst, emChromaticity, exEst, exCoeffs, predRefl, predFl, hist ] = fiRecReflAndFlMultistep( measVals, cameraMat, cameraGain, cameraOffset, illuminant, basisRefl, DB, basisEx, alpha, gamma, varargin )

% [ reflEst, reflCoeffs, emEst, emChromaticity, exEst, exCoeffs, predRefl, predFl, hist ] = fiRecReflAndFlMultistep( measVals, cameraMat, cameraGain, cameraOffset, illuminant, basisRefl, DB, basisEx, alpha, gamma, ... )
%
% Perform reflectance and fluorescence excitation, emission spectra
% estimation using the multistep approach of Fu et al. decsribed in
% 'Reflectance and Fluorescent Spectra Recovery based on Fluorescent
% Chromaticity Invariance under Varying Illumination,' CVPR 2014.
%
% Inputs:
%    measVals - a (f x c x s) matrix containing pixel intensities of s 
%      different surfaces captured with f camera channels captured under c
%      different illuminants.
%    cameraMat - a (w x c) matrix containing the spectral responsivity
%      functions of the c camera channels sampled at w wavebands.
%    cameraGain - a (f x c x s) matrix of linear camera model gains for each
%      filter-channel-surface combination.
%    cameraOffset - a (f x c x s) matrix of linear camera model offsets for
%      each filter-channel-surface combination.
%    illuminant - a (w x c) matrix containing the spectral power
%      distributions of the c illuminants. 
%    basisRefl, basisEx - (w x n) matrices of n linear basis
%      functions representing reflectance, and excitation spectra
%      respectively.
%    DB - a (w x t) array containing a collection of t emission spectra. The 
%      emission spectrum estimate will be selected from this data set.
%    alpha - scalar tuning parameter controlling the smoothness of
%      reflectance estimates.
%    gamma - scalar tuning parameter controlling the smoothness of
%      fluorescence excitation estimates.
%
% Inputs (optional):
%    'maxIter' - a scalar defining the maximal number of biconvex solver
%       iterations (default = 100).
%    'eps' - biconvex solver convergence threshold (default = 1e-8).
%
% Outputs:
%    reflEst - a (w x s) matrix of estimated surface spectral reflectances.
%    rfCoeffs - a (n x s) matrix expressing the estimated surface spectral 
%      reflectances in terms of the linear basis weights.
%    emEst - a (w x s) matrix of the estimated surface emission spectra.
%    emChromaticity - a (f x s) matrix expressing the estimated surface emission 
%      spectra chromaticities.
%    exEst - a (w x s) matrix of the estimated surface excitation spectra.
%    exCoeffs - a (n x s) matrix expressing the estimated surface excitation 
%      spectra in terms of the linear basis weights.
%    predRefl - a (f x c x s) matrix of values representing the reflected
%      light intensities for each filter-illuminant-surface combination.
%    predFl - a (f x c x s) matrix of values representing the fluoresced
%      light intensities for each filter-illuminant-surface combination.
%      Specifically, for ideal, noiseless measurements the following holds
%              predRefl + predFl = (measVals - cameraOffset)
%    hist - a s-dimensional cell array containing the objective function
%      values at successive minimization steps of the biconvex solver.
%
% Copyright, Henryk Blasinski 2016.


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
    
    % First estimate reflectance and fluorescence emission chromaticity. 
    % See section 3.2 in the original paper
    [reflEst(:,i), reflCoeffs(:,i), emChromaticity(:,i), hist{i}] = fiRecOneRefl(input,...
        cameraMat,cameraGain(:,:,i),illuminant,basisRefl,alpha, varargin{:});
    
    % Estimate the excitation spectrum (section 3.3)
    [exEst(:,i), exCoeffs(:,i)] = fiRecOneEx(input,reflCoeffs(:,i),emChromaticity(:,i),...
        cameraMat,cameraGain(:,:,i),illuminant,basisRefl,basisEx,gamma);
    
    % Estimate the emission spectrum (section 3.4)
    [emEst(:,i)] = fiRecOneEm(emChromaticity(:,i),cameraMat,DB);
    
    predRefl(:,:,i) = cameraGain(:,:,i).*(cameraMat'*diag(reflEst(:,i))*illuminant);
    
    flEx = exEst(:,i)'*illuminant;
    flEm = cameraMat'*emEst(:,i);
    
    % The gain parameter estimation is not described in the paper.
    
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
    
end


end