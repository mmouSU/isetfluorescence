function [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, predRefl, predFl, hist  ] = fiRecOneReflAndFl( measVals, cameraMat, cameraGain, illuminant, basisRefl, basisEm, basisEx, alpha, beta, gamma, varargin )


p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('forceIter',false,@islogical);
p.addParamValue('epsilon',1e-8);
p.addParamValue('maxIter',100);
p.parse(varargin{:});
inputs = p.Results;


nEmBasis = size(basisEm,2);
nExBasis = size(basisEx,2);
nReflBasis = size(basisRefl,2);
nWaves = size(illuminant,1);
nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);

% Create the camera linear model
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


% Begin a biconvex problem iteration
hist.objValsReEm = zeros(inputs.maxIter+1,1);
hist.objValsReEx = zeros(inputs.maxIter+1,1);
for i=1:inputs.maxIter
    
     
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
        
    % Exit if the objective improvement is small.
    if abs(hist.objValsReEm(i) - hist.objValsReEx(i)) <= inputs.epsilon
        break;
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
% max(ex)=1, and all the scaling is incorporated in the emission.



DM = tril(emEst*exEst',-1);
[U, S, V] = svd(DM);

emEst = U(:,1)*sign(min(U(:,1)));
exEst = S(1,1)*V(:,1)*sign(min(V(:,1)));

sf = max(exEst);
exEst = exEst/sf;
emEst = emEst*sf;

% Compute the prediction
predRefl = cameraGain.*(cameraMat*diag(reflEst)*illuminant);
predFl = cameraGain.*(cameraMat*DM*illuminant);


end

