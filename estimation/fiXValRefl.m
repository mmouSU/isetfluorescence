function [ optLambda, error, reflEstimate, reflCoeffs, predResp ] = fiXValRefl( measVals, cameraResp, cameraGain, cameraOffset, illuminant, basisFcns, lambdaSet, reflReference )

% [ optLambda, convMat, reflEstimate, reflCoeffs, predResp ] = fiXValRefl( measVals, cameraResp, cameraGain, cameraOffset, illuminant, basisFcns, lambdaSet, reflReference )
%
% This function is used to perform the cross-validation search for the
% optimal smoothness parameter lambda in the fiRecRefl
% function.
%
% Inputs:
%   see help for fiRecRefl
%   lambdaSet - a vector containing the values of lambda used in the cross
%   validation process
%   reflReference - a nWaves x nSamples matrix of reflectance reference
%
% Outputs:
%   optLambda - a the optimal value of the smoothness parameter
%   error - a vector of mean RMS error values for each value of the
%   smoothness parameter from the lambdaSet
%   see help for linearFindReflectanceV2
%
% Copyright, Henryk Blasinski 2014

nLambda = length(lambdaSet);

measVals = measVals - cameraOffset;

error = Inf*ones(nLambda+1,1);
optLambda = [];
reflEstimate = [];
predResp = [];
reflCoeffs = [];

nrm = @(x) sqrt(mean((x).^2));

for l=1:nLambda
   
    fprintf('Lambda: %5.5f (%5i/%i); rms = ',lambdaSet(l),l,nLambda);
    
    
    [est, coeffs, pred] = fiRecRefl(measVals,cameraResp,cameraGain,illuminant,basisFcns,lambdaSet(l));

    
    error(l+1) = mean(nrm(est - reflReference));
    if error(l+1) <=  min(error(1:l-1))
        optLambda = lambdaSet(l);
        reflEstimate = est;
        predResp = pred;
        reflCoeffs = coeffs;
    end
    
    fprintf(' %f\n',error(l));
end
error = error(2:end);

end

