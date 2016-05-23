function [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, predRefl, predFl, hist  ] = fiRecReflAndFl( measVals, camera, cameraGain, cameraOffset, illuminant, basisRefl, basisEm, basisEx, alpha, beta, gamma, varargin )

% [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, predRefl, predFl, hist  ] = fiRecReflAndFl( measVals, camera, cameraGain, cameraOffset, illuminant, basisRefl, basisEm, basisEx, alpha, beta, gamma, ... )
%
% This is a wrapper function for single fluorophore estimation algorithm to
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
%    basisRefl, basisEx, basisEm - (w x n) matrices of n linear basis
%      functions representing reflectance, excitation and emission spectra
%      respectively
%    alpha - scalar tuning parameter controlling the smoothness of
%      reflectance estimates
%    beta - scalar tuning parameter controlling the smoothness of
%      fluorescence excitation estimates
%    gamma - scalar tuning parameter controlling the smoothness of
%      fluorescence emission estimates
%
% Inputs (optional):
%    'reflRef' - a (w x s) matrix of reference surface reflectances. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure. 
%    'exRef' - a (w x s) matrix of reference fluorophore excitation spectra. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure. 
%    'exRef' - a (w x s) matrix of reference fluorophore emission spectra. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure.
%    'pixelRef' - a boolean value indicating if the error between predicted
%      and mesured pixel intensities is to be computed at every iteration.
%      (default = false).
%    Additioanl inputs control the bi-convex optimizer are transferred directly
%    to the estimation function directly. See fiRecOneReflAndFl for details.
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


p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('reflRef',[]);
p.addParamValue('exRef',[]);
p.addParamValue('emRef',[]);
p.parse(varargin{:});
inputs = p.Results;

nSamples = size(measVals,3);
nFilters = size(camera,2);
nChannels = size(illuminant,2);
nWaves = size(camera,1);
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
    
    if ~isempty(inputs.exRef)
        exRef = inputs.exRef(:,i);
    else
        exRef = [];
    end
    if ~isempty(inputs.emRef)
        emRef = inputs.emRef(:,i);
    else
        emRef = [];
    end
    
    if ~isempty(inputs.reflRef);
        reflRef = inputs.reflRef(:,i);
    else
        reflRef = [];
    end
    
    [reflEst(:,i), rfCoeffs(:,i), emEst(:,i), emCoeffs(:,i), exEst(:,i), exCoeffs(:,i), predRefl(:,:,i), predFl(:,:,i), hist{i}] = ...
        fiRecOneReflAndFl(input,camera,cameraGain(:,:,i),illuminant,basisRefl,basisEm,basisEx,alpha,beta,gamma,varargin{:},'reflRef',reflRef,...
        'exRef',exRef,'emRef',emRef);

    nIter = length(hist{i}.objValsReEm);

    fprintf('Done (%i iterations)\n',nIter);

end


end

