function [ reflEst, emEst, exEst, dMatEst, predRefl, predFl, hist  ] = fiRecReflAndFlNucNorm( measVals, cameraMat, cameraGain, cameraOffset, illuminant, alpha, sigma, varargin )

% [ reflEst, emEst, exEst, dMatEst, predRefl, predFl, hist  ] = fiRecOneReflAndFlNucNorm( measVals, cameraMat, cameraGain, cameraOffset, illuminant, alpha, sigma, ... )
%
% Estimate the reflectance and fluorescence of a collection of surfaces
% using an adaptation of the algorithm from Suo et al. 'Bispectral coding:
% compressive and high0quality acquisition of fluorescence and
% reflectance,' Optics Express 2014. In this implementation we solve for
% spectral properties (F and R) of a single pixel, i.e. we solve the
% optimization problem (5) using our own ADMM solver. Our attemnt was to
% follow this paper to the letter. Specifically we do not impose any
% additional constratins that could follow from the physics of the problem.
% For example the matrix R is not restricted to being diagonal nor
% non-negative.
%
% Inputs:
%   measVals - a (f x c x s) array of pixel measurements through f filters,
%      under c illuminants and for s surfaces
%   cameraMat - a (w x f) array of spectral responsivity
%      functions of the f camera filters.
%   cameraGain - a (f x c x s) array of camera gains applied to
%      every channel-filter-surface combination.
%   cameraOffset - a (f x c x s) array of camera offsets applied to
%      every channel-filter-surface combination.
%   illuminant - a (w x c) array of the illuminant spectral
%      power distributions of the c illuminants.
%   alpha - a scalar controling the l1 norm applied to the reflectance
%      matrix R.
%   sigma - a scalar controlling the admissible error in the predicted
%      pixel values.
%
% Inputs (optional):
%   'maxIter' - the maximal number of ADMM interations (default = 500).
%   'rescaleRho' - implement a heuristic algorithm improving ADMM
%      convergence from Boyd 2011 (default = true).
%   'tauIncr' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 2).
%   'tauDecr' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 2).
%   'rhoInit' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 1).
%   'mu' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 10).
%   'verbose' - display residual values every 10th iteration 
%      (default = true).
%   'epsilon' - the convergence threshold for the ADMM algorithm 
%      (default = 1e-6).
%
% Outputs:
%   reflEst - an (s x 1) cell array of square matrices representing the
%      surface spectral reflectances. Note the algorithm does not impose 
%      any constraints on the elements of reflEst. The reflectance estimates
%      do not have to be diagonal or non-negative.
%   emEst - an (w x s) array of emission spectra estimates for the s surfaces.
%   exEst - an (w x s)  array of excitation spectra estimates for the s surfaces.
%   F - an (s x 1) cell array of the Donaldson matrix estimates.
%   predRefl - a (f x c x s) array of predicted pixel
%      intensities due to reflected light component.
%   predFl - a (f x c x s) array of predicted pixel intensities due to 
%      fluoresced light component.
%   hist - a (s x 1) array of structures containing primal and dual residuals
%      of the ADMM solver, as well as number of iterations and errors of the PCG
%      algorithm used to solve least-squares minimization at every ADMM
%      iteration.
% 
% Copyright, Henryk Blasinski 2016.


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

    fprintf('Processing sample %i ... ',i);
    
    input = measVals(:,:,i) - cameraOffset(:,:,i);
    
    [reflEst{i}, emEst(:,i), exEst(:,i), dMatEst{i}, predRefl(:,:,i), predFl(:,:,i), hist{i}] = ...
    fiRecOneReflAndFlNucNorm(input,cameraMat,cameraGain(:,:,i),illuminant,alpha,sigma,varargin{:});
    
    fprintf('Done!\n');
end

end

