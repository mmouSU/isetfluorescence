function [ reflEst, rfCoeffs, emEst, emCoeffs, emWghts, predRefl, predFl, hist  ] = fiRecReflAndEm( measVals, cameraMat, cameraGain, cameraOffset, illuminant, basisRefl, basisEm, alpha, beta, varargin )

% [ reflEst, rfCoeffs, emEst, emCoeffs, emWghts, predRefl, predFl, hist  ] = fiRecReflAndEm( measVals, camera, cameraGain, cameraOffset, illuminant, basisRefl, basisEm, basisEx, alpha, beta, ... )
%
% This is a wrapper function for the CIM estimation algorithm to
% perform estimation on multiple surfaces measured with the same imaging system
% Measurements of different surfaces share the imaging system properties and 
% the tuning parameter settings, but are otherwise independent (no spatial 
% smoothness constraints). 
%
% Inputs (required):
%    measVals - a (f x c x s) matrix containing pixel intensities of s 
%      different surfaces captured with f camera channels captured under c
%      different illuminants.
%    camera - a (w x c) matrix containing the spectral responsivity
%      functions of the c camera channels sampled at w wavebands.
%    cameraGain - a (f x c x s) matrix of linear camera model gains for each
%      filter-channel-surface combination.
%    cameraOffset - a (f x c x s) matrix of linear camera model offsets for
%      each filter-channel-surface combination
%    illuminant - a (w x c) matrix containing the spectral power
%      distributions of the c illuminants. 
%    basisRefl, basisEm - (w x n) matrices of n linear basis
%      functions representing reflectance, excitation and emission spectra
%      respectively
%    alpha - scalar tuning parameter controlling the smoothness of
%      reflectance estimates
%    beta - scalar tuning parameter controlling the smoothness of
%      fluorescence excitation estimates
%
% Inputs (optional):
%    'tol' - a scalar describing the change in the objective function
%      value that causes the iterative algorithm to terminate (default: 1e-8)
%    'maxIter' - maximal number of iterations of the biconvex algorithm. Two
%      convex problems are solved for every iteration: one to estimate
%      reflectance and emission spectrum, the other to estimate reflectance
%      and emission scale (default 100).
%
% Outputs:
%    reflEst - a (w x s) matrix of estimated surface spectral reflectances.
%    rfCoeffs - a (n x s) matrix expressing the estimated surface spectral 
%      reflectances in terms of the linear basis weights.
%    emEst - a (w x s) matrix of the estimated surface emission spectra.
%    emCoeffs - a (n x s) matrix expressing the estimated surface emission 
%      spectra in terms of the linear basis weights.
%    emWghts - a (c x s) array of scale factors representing for every surface
%      the intensity of fluorescence emission under each illuminant.
%    predRefl - a (f x c x s) matrix of values representing the reflected
%      light intensities for each filter-illuminant-surface combination.
%    predFl - a (f x c x s) matrix of values representing the fluoresced
%      light intensities for each filter-illuminant-surface combination.
%      Specifically, for ideal, noiseless measurements the following holds
%              predRefl + predFl = (measVals - cameraOffset)
%    hist - a s-dimensional cell array containing the objective function
%      values at successive minimization steps.
%
% Copyright, Henryk Blasinski 2016.

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

