function [ reflEst, rfCoeffs, emEst, emCoeffs, emWghts, predRefl, predFl, hist  ] = fiRecOneReflAndEm( measVals, cameraMat, cameraGain, illuminant, basisRefl, basisEm, alpha, beta, varargin )

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

for i=1:inputs.maxIter;
    
    
    
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
    
    % Exit if the objective improvement is small.
    if abs(hist.objValsReWe(i) - hist.objValsReEm(i)) <= inputs.tol,
        break;
    end
    
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

