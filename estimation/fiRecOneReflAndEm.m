function [ reflEst, rfCoeffs, emEst, emCoeffs, emWghts, predRefl, predFl, hist  ] = fiRecOneReflAndEm( measVals, cameraMat, cameraGain, illuminant, basisRefl, basisEm, alpha, beta, varargin )

% [ reflEst, rfCoeffs, emEst, emCoeffs, emWghts, predRefl, predFl, hist  ] = fiRecOneReflAndEm( measVals, cameraMat, cameraGain, illuminant, basisRefl, basisEm, alpha, beta, ... )

% This function implements the color invariant (CIM) fluorophore model and 
% estimates the reflectance and fluorescence emission spectra 
% from pixel intensities using a bi-convex estimation algorithm. This 
% implementation contains tuning parameters to independently adjust the 
% smoothness of the excitation and emission spectra. 
%
% Inputs (required):
%    measVals - a (f x c) matrix containing pixel intensities of a 
%      surface captured with f camera channels and under c
%      different illuminants.
%    camera - a (w x c) matrix containing the spectral responsivity
%      functions of the c camera channels sampled at w wavebands.
%    cameraGain - a (f x c) matrix of linear camera model gains for each
%      filter-channel combination.
%    basisRefl, basisEm - a (w x n) matrices of c linear basis
%      functions representing reflectance and emission spectra
%      respectively
%    alpha - scalar tuning parameters controlling the smoothness of
%      reflectance estimates
%    beta - scalar tuning parameters controlling the smoothness of
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
%    reflEst - a (w x 1) vector of the estimated surface spectral reflectance.
%    rfCoeffs - a (s x 1) vector expressing the estimated surface spectral 
%      reflectance in terms of the linear basis weights.
%    emEst - a (w x 1) vector of the estimated emission spectrum.
%    emCoeffs - a (s x 1) vector expressing the estimated surface emission 
%      spectrum in terms of the linear basis weights.
%    emWghts - a (c x 1) vector of scale factors representing the intensity
%      of fluorescence emission under each illuminant.
%    predRefl - a (f x c) matrix of values representing the reflected
%      light intensities for each filter-illuminant combination.
%    predFl - a (f x c) matrix of values representing the fluoresced
%      light intensities for each filter-illuminant combination.
%      Specifically, for ideal, noiseless measurements the following holds
%              predRefl + predFl = (measVals - cameraOffset)
%    hist - structure containing the objective function values at 
%      successive minimization steps.
%
% Copytight, Henryk Blasinski 2016




p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('tol',1e-6);
p.addParamValue('maxIter',100);
p.parse(varargin{:});
inputs = p.Results;


% This function assumes a fluorescent model where the emission spectrum is
% fixed and only it's intensity varies.
% See Separating Reflective and Fluorescent Components Using High Frequency
% Illumination in the Spectral Domain.

nFlBasis = size(basisEm,2);
nReflBasis = size(basisRefl,2);
nWaves = size(illuminant,1);
nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);

scaleFactor = max(cameraGain(:));
cameraGain = cameraGain/scaleFactor;
illuminant = illuminant*scaleFactor;
cameraMat = cameraMat';


% Create a penalty on roughness
R = [diag(ones(nWaves-1,1)) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) diag(ones(nWaves-1,1))];
Rem = R*basisEm;
Rrefl = R*basisRefl;


% Assume some initial distribution of variables. 
emWghts = ones(nChannels,1);
rfCoeffs = ones(nReflBasis,1);
emCoeffs = ones(nFlBasis,1);

% Begin a biconvex problem iteration
hist.objValsReEm = zeros(inputs.maxIter,1);
hist.objValsReWe = zeros(inputs.maxIter,1);

for i=1:inputs.maxIter
    
    tic;
    
    % Assume that you know the fluorescence weights, search for the shape
    % of the fluorescence and reflectance
    cvx_begin quiet
        cvx_precision default
        variables rfCoeffs(nReflBasis,1) emCoeffs(nFlBasis,1)
        minimize sum(sum_square_abs(measVals - cameraGain.*(cameraMat*diag(basisRefl*rfCoeffs)*illuminant) - cameraGain.*(cameraMat*(basisEm*emCoeffs*(emWghts'))))) + alpha*norm(Rrefl*rfCoeffs,2) + beta*norm(Rem*emCoeffs,2)
        subject to
            1 >= basisRefl*rfCoeffs >= 0
            basisEm*emCoeffs >= 0
    cvx_end
    if strcmp(cvx_status,'Failed') == 1, break; end
    
    hist.objValsReEm(i) = cvx_optval;
        
    %% Re-arrange the fluorescence component
    
    cvx_begin quiet
        cvx_precision default
        variables emWghts(nChannels,1) rfCoeffs(nReflBasis,1)
        minimize sum(sum_square_abs(measVals - cameraGain.*(cameraMat*diag(basisRefl*rfCoeffs)*illuminant) - cameraGain.*(cameraMat*(basisEm*emCoeffs*(emWghts'))))) + alpha*norm(Rrefl*rfCoeffs,2) + beta*norm(Rem*emCoeffs,2)
        subject to
            emWghts >= 0
            1 >= basisRefl*rfCoeffs >= 0
    cvx_end
    if strcmp(cvx_status,'Failed') == 1, break; end

    
    hist.objValsReWe(i) = cvx_optval;
    
    tend = toc;
    % Exit if the objective improvement is small.
    if abs(hist.objValsReWe(i) - hist.objValsReEm(i)) <= inputs.tol
        break;
    end
    
    hist.computeTime(i) = tend;
    
end

hist.objValsReWe = hist.objValsReWe(1:i);
hist.objValsReEm = hist.objValsReEm(1:i);

% Estimate Fluorescence shape.
emEst = basisEm*emCoeffs;

scaleFac = sum(emEst);
emEst = emEst/scaleFac;

emWghts = emWghts*scaleFac;

% Estimate Reflectnace
reflEst = basisRefl*rfCoeffs;

% Compute the prediction
predRefl = cameraGain.*(cameraMat*(diag(reflEst)*illuminant));
predFl = cameraGain.*(cameraMat*(emEst*emWghts'));



end

