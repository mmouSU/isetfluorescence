function [ reflEst, reflCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMat, predRefl, predFl, hist  ] = fiRecReflAndMultiFl( measVals, cameraMat, illuminant, cameraGain, cameraOffset, basisRefl, basisEm, basisEx, alpha, beta, gamma, eta, varargin)

% [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMat, predRefl, predFl, hist  ] = fiRecReflAndMultiFl( measVals, cameraMat, illuminant, cameraGain, cameraOffset, basisRefl, basisEm, basisEx, alpha, beta, gamma, eta, varargin)
%
% This is a wrapper function for the multi-fluorophore estimation algorithm to
% perform estimation on multiple surfaces measured with the same imaging system
% Measurements of different surfaces share the imaging system properties and 
% the tuning parameter settings, but are otherwise independent (no spatial 
% smoothness constraints).
%
% Typically beta=gamma, i.e. the smothness penalty on excitation and
% emission are equal.
%
% Inputs (required):
%    measVals - a (f x c x s) matrix containing pixel intensities of s 
%      different surfaces captured with f camera channels captured under c
%      different illuminants.
%    cameraMat - a (w x c) matrix containing the spectral responsivity
%      functions of the c camera channels sampled at w wavebands.
%    illuminant - a (w x s) matrix containing the spectral power
%      distributions of s illuminant energies 
%    cameraGain - a (f x c x s) matrix of linear camera model gains for each
%      filter-channel-surface combination.
%    cameraOffset - a (f x c x s) matrix of linear camera model offsets for
%      each filter-channel-surface combination
%    basisRefl, basisEx, basisEm - a (w x n) matrices of c linear basis
%      functions representing reflectance, excitation and emission spectra
%      respectively
%    alpha - scalar tuning parameter controlling the smoothness of
%      reflectance estimates
%    beta - scalar tuning parameter controlling the smoothness of
%      fluorescence excitation estimates
%    gamma - scalar tuning parameter controlling the smoothness of
%      fluorescence emission estimates
%    eta -  a scalar tuning parameter controlling the nuclear norm (i.e.
%      the number of fluorophores) of the Donaldson matrixestimate.
%
% Inputs (optional):
%    'reflRef' - a (w x s) matrix of reference surface reflectances. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure. 
%    'dMatRef' - a (s x 1) cell array of reference Donaldson matrices. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure. 
%    'pixelRef' - a boolean value indicating if the error between predicted
%      and mesured pixel intensities is to be computed at every iteration.
%      (default = false).
%    Additioanl inputs control the ADMM optimizer are transferred directly
%    to the estimation function directly. See fiRecOneReflAndMultiFl for 
%    details.
%
% Outputs:
%    reflEst - a (w x s) matrix of estimated surface spectral reflectances.
%    rfCoeffs - a (n x s) matrix expressing the estimated surface spectral 
%      reflectances in terms of the linear basis weights.
%    emEst - a (w x s) matrix of the estimated surface emission spectra.
%    emCoeffs - a (n x s) matrix expressing the estimated surface emission 
%      spectra in terms of the linear basis weights.
%    exEst - a (w x s) matrix of the estimated surface excitation spectra.
%    exCoeffs - a (n x s) matrix expressing the estimated surface excitation 
%      spectra in terms of the linear basis weights.
%    dMat - a (s x 1) cell array of Donaldson matrix estimates.
%    predRefl - a (f x c x s) matrix of values representing the reflected
%      light intensities for each filter-illuminant-surface combination.
%    predFl - a (f x c x s) matrix of values representing the fluoresced
%      light intensities for each filter-illuminant-surface combination.
%      Specifically, for ideal, noiseless measurements the following holds
%              predRefl + predFl = (measVals - cameraOffset)
%    hist - a c-dimensional cell array containing the objective function
%      values at successive minimization steps.
%
% Copyright, Henryk Blasinski 2016.


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

reflEst = zeros(nWaves,nSamples);
reflCoeffs = zeros(nBasisRefl,nSamples);

emEst = zeros(nWaves,nSamples);
emCoeffs = zeros(nBasisEm,nSamples);

exEst = zeros(nWaves,nSamples);
exCoeffs = zeros(nBasisEx,nSamples);

dMat = cell(nSamples,1);

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

    [reflEst(:,i), reflCoeffs(:,i), emEst(:,i), emCoeffs(:,i), exEst(:,i), exCoeffs(:,i), dMat{i}, predRefl(:,:,i), predFl(:,:,i), hist{i}] = ...
    fiRecOneReflAndMultiFl(input,cameraMat,illuminant,cameraGain(:,:,i),basisRefl,basisEm,basisEx,alpha,beta,gamma,eta,varargin{:},'dMatRef',dMatRef,'reflRef',reflRef);
    
    fprintf('Done!\n');
end

end

