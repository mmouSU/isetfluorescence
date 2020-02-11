function [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, predRefl, predFl, hist  ] = fiRecOneReflAndFl( measVals, cameraMat, cameraGain, illuminant, basisRefl, basisEm, basisEx, alpha, beta, gamma, varargin )

% [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, predRefl, predFl, hist  ] = fiRecOneReflAndFl( measVals, cameraMat, cameraGain, illuminant, basisRefl, basisEm, basisEx, alpha, beta, gamma, ... )

% This function implements the single fluorophore model and estimates the
% reflectance, fluorescence excitation and emission spectra from pixel 
% intensities using a bi-convex estimation algorithm. This implementation
% contains tuning parameters to independently adjust the smoothness of the
% excitation and emission spectra. 
%
% Inputs (required):
%    measVals - a (f x c) matrix containing pixel intensities of a 
%      surface captured with f camera channels and under c
%      different illuminants.
%    camera - a (w x c) matrix containing the spectral responsivity
%      functions of the c camera channels sampled at w wavebands.
%    cameraGain - a (f x c) matrix of linear camera model gains for each
%      filter-channel combination.
%    basisRefl, basisEx, basisEm - a (w x n) matrices of c linear basis
%      functions representing reflectance, excitation and emission spectra
%      respectively
%    alpha - scalar tuning parameters controlling the smoothness of
%      reflectance estimates
%    beta - scalar tuning parameters controlling the smoothness of
%      fluorescence excitation estimates
%    gamma - scalar tuning parameters controlling the smoothness of
%      fluorescence emission estimates
%
% Inputs (optional):
%    'epsilon' - a scalar describing the change in the objective function
%      value that causes the iterative algorithm to terminate (default: 1e-8)
%    'maxIter' - maximal number of iterations of the biconvex algorithm. Two
%      convex problems are solved for every iteration: one to estimate
%      reflectance and emission spectrum, the other to estimate reflectance
%      and excitation spectrum (default 100).
%    'reflRef' - a (w x 1) vector of reference surface reflectance. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure. 
%    'exRef' - a (w x 1) vector of reference fluorophore excitation spectrum. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure. 
%    'exRef' - a (w x 1) matrix of reference fluorophore emission spectrum. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure.
%    'pixelRef' - a boolean value indicating if the error between predicted
%      and mesured pixel intensities is to be computed at every iteration.
%      (default = false).
%
% Outputs:
%    reflEst - a (w x 1) vector of the estimated surface spectral reflectance.
%    rfCoeffs - a (s x 1) vector expressing the estimated surface spectral 
%      reflectance in terms of the linear basis weights.
%    emEst - a (w x 1) vector of the estimated emission spectrum.
%    emCoeffs - a (s x 1) vector expressing the estimated surface emission 
%      spectrum in terms of the linear basis weights.
%    exEst - a (w x 1) vector of the estimated surface excitation spectrum.
%    exCoeffs - a (s x 1) vector expressing the estimated surface excitation 
%      spectrum in terms of the linear basis weights.
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
p.addParamValue('epsilon',1e-8);
p.addParamValue('maxIter',100);
p.addParamValue('reflRef',[]);
p.addParamValue('exRef',[]);
p.addParamValue('emRef',[]);
p.addParamValue('pixelRef',false);
p.parse(varargin{:});
inputs = p.Results;


nEmBasis = size(basisEm,2);
nExBasis = size(basisEx,2);
nReflBasis = size(basisRefl,2);
nWaves = size(illuminant,1);
nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);

% Scale the inputs to improve numerical stability of the solution.
cameraMat = cameraMat';
scaleFac = max(cameraGain(:));
cameraGain = cameraGain/scaleFac;
illuminant = illuminant*scaleFac;


% Create a penalty on roughness
R = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];
Rem = R*basisEm;
Rrefl = R*basisRefl;
Rex = R*basisEx;

% Assume some initial distribution of variables. 
exCoeffs = ones(nExBasis,1)/nExBasis;


hist.objValsReEm = zeros(inputs.maxIter+1,1);
hist.objValsReEx = zeros(inputs.maxIter+1,1);

if inputs.pixelRef
    hist.pixelErr = zeros(inputs.maxIter+1,1);
end
if ~isempty(inputs.reflRef)
    hist.reflErr = zeros(inputs.maxIter+1,1);
end
if ~isempty(inputs.exRef)
    hist.exErr = zeros(inputs.maxIter+1,1);
end
if ~isempty(inputs.emRef)
    hist.emErr = zeros(inputs.maxIter+1,1);
end


% Begin a biconvex problem iteration
for i=1:inputs.maxIter
    
    tic;
     
    cvx_begin quiet
        cvx_precision default
        variables rfCoeffs(nReflBasis,1) emCoeffs(nEmBasis,1)
        minimize sum(sum_square_abs(measVals - cameraGain.*(cameraMat*diag(basisRefl*rfCoeffs)*illuminant) - cameraGain.*(cameraMat*tril((basisEm*emCoeffs)*(basisEx*exCoeffs)',-1)*illuminant))) + alpha*norm(Rrefl*rfCoeffs,2) + beta*norm(Rem*emCoeffs,2) + gamma*norm(Rex*exCoeffs,2)
        subject to
            1 >= basisRefl*rfCoeffs >= 0
            basisEm*emCoeffs >= 0
    cvx_end
    if strcmp(cvx_status,'Failed') == 1, break; end
    
    hist.objValsReEm(i) = cvx_optval;
    
    cvx_begin quiet
        cvx_precision default
        variables rfCoeffs(nReflBasis,1) exCoeffs(nExBasis,1)
        minimize sum(sum_square_abs(measVals - cameraGain.*(cameraMat*diag(basisRefl*rfCoeffs)*illuminant) - cameraGain.*(cameraMat*tril((basisEm*emCoeffs)*(basisEx*exCoeffs)',-1)*illuminant))) + alpha*norm(Rrefl*rfCoeffs,2) + beta*norm(Rem*emCoeffs,2) + gamma*norm(Rex*exCoeffs,2)
        subject to
            1 >= basisRefl*rfCoeffs >= 0
            basisEx*exCoeffs >= 0
    cvx_end
    if strcmp(cvx_status,'Failed') == 1, break; end
   
    hist.objValsReEx(i) = cvx_optval;
    hist.computeTime(i) = toc;
    
    % Exit if the objective improvement is small.
    if abs(hist.objValsReEm(i) - hist.objValsReEx(i)) <= inputs.epsilon
        break;
    end
    
    %% Copute errors (if the reference is provided), 
    % This computation this slows the algorithm a bit.
    
    % Estimate Reflectnace
    if ~isempty(inputs.reflRef)
        reflEst = basisRefl*(rfCoeffs);
        hist.reflErr(i) = fiComputeError(reflEst,inputs.reflRef,'absolute');
    end

    % Estimate Fluorescence 
    if ~isempty(inputs.exRef)
        exEst = basisEx*exCoeffs;
        hist.exErr(i) = fiComputeError(exEst,inputs.exRef,'normalized');
    end
    if ~isempty(inputs.emRef)
        emEst = basisEm*emCoeffs;
        hist.emErr(i) = fiComputeError(emEst,inputs.emRef,'normalized');
    end

    % Estimate pixels
    if inputs.pixelRef
        emEst = basisEm*emCoeffs;
        exEst = basisEx*exCoeffs;

        reflEst = basisRefl*(rfCoeffs);
        
        DM = tril(emEst*exEst',-1);
        % Compute pixel intensities predictions
        predRefl = cameraGain.*(cameraMat*diag(reflEst)*illuminant);
        predFl = cameraGain.*(cameraMat*DM*illuminant);
        
        hist.pixelErr(i) = fiComputeError(predRefl(:) + predFl(:),measVals(:),'');
    end
    
   
end

hist.objValsReEm = hist.objValsReEm(1:i);
hist.objValsReEx = hist.objValsReEx(1:i);

% Estimate Reflectnace
reflEst = basisRefl*(rfCoeffs);

% Estimate Fluorescence shape
emEst = basisEm*emCoeffs;
exEst = basisEx*exCoeffs;

% Now we have to re-scale the data. We follow the convention that
% max(ex)=1, and all the scaling is incorporated in the emission. Also note
% that we use SVD to estimate the excitation and emission spectra.
% Technically we could return the values given directly by the weights, 
% but since we explicitly impose lower triangularity, these results 
% sometimes include 'bumps' that violate this condition. 

DM = tril(emEst*exEst',-1);
[U, S, V] = svd(DM);

emEst = U(:,1)*sign(min(U(:,1)));
exEst = S(1,1)*V(:,1)*sign(min(V(:,1)));

sf = max(exEst);
exEst = exEst/sf;
emEst = emEst*sf;

% Compute pixel intensities predictions
predRefl = cameraGain.*(cameraMat*diag(reflEst)*illuminant);
predFl = cameraGain.*(cameraMat*DM*illuminant);


end

