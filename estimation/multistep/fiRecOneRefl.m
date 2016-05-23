function [ reflEst, reflCoeffs, E, hist ] = fiRecOneRefl( measVals, cameraMat, cameraGain, illuminant, basisRefl, alpha, varargin )

% [ reflEst, reflCoeffs, E, hist ] = fiRecOneRefl( measVals, cameraMat, cameraGain, illuminant, basisRefl, alpha, varargin )
%
% Estimate reflectance and fluorescence emission chromaticity. This function 
% implements the approach from the Fu et al. CVPR 2014 paper, section 3.2
%
% Inputs:
%    measVals - a (f x c x s) matrix containing pixel intensities of s 
%      different surfaces captured with f camera channels captured under c
%      different illuminants.
%    cameraMat - a (w x c) matrix containing the spectral responsivity
%      functions of the c camera channels sampled at w wavebands.
%    cameraGain - a (f x c x s) matrix of linear camera model gains for each
%      filter-channel-surface combination.
%    illuminant - a (w x c) matrix containing the spectral power
%      distributions of the c illuminants. 
%    basisRefl, basisEx - (w x n) matrices of n linear basis
%      functions representing reflectance, and excitation spectra
%      respectively.
%    alpha - scalar tuning parameter controlling the smoothness of
%      reflectance estimates.
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
%    E - a (f x s) array of fluorescence emission chromaticities. 
%    hist - a s-dimensional cell array containing the objective function
%      values at successive minimization steps of the biconvex solver.
%
% Copyright, Henryk Blasinski 2016.

p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('maxIter',100);
p.addParamValue('eps',1e-8);
p.parse(varargin{:});

inputs = p.Results;

nReflBasis = size(basisRefl,2);
nWaves = size(illuminant,1);
nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);

q_m_nj = zeros(nFilters,nChannels,nReflBasis);
for f=1:nFilters
    for c=1:nChannels
        q_m_nj(f,c,:) = cameraGain(f,c)*(cameraMat(:,f)'*diag(illuminant(:,c))*basisRefl);
    end
end

R = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1) eye(nWaves-1,1)];

E = mean(measVals,2);
E = E/sum(E);

hist.objValsRefl = zeros(inputs.maxIter,1);
hist.objValsChr = zeros(inputs.maxIter,1);

for i=1:inputs.maxIter
    
    % Estimate reflectance
    cvx_begin quiet
        variable reflCoeffs(nReflBasis,1)
        coeffQ = sum(q_m_nj.*repmat(shiftdim(reflCoeffs,-2),[nFilters,nChannels,1]),3);
        minimize sum(sum_square(measVals - E*sum(measVals - coeffQ) - coeffQ)) +  alpha*sum(sum_square(R*basisRefl*reflCoeffs))
        subject to
            basisRefl*reflCoeffs >= 0
    cvx_end
    
    hist.objValsRefl(i) = cvx_optval;
    
    % Estimate fluorescence emission chromaticity 
    cvx_begin quiet
        variable E(nFilters,1)
        coeffQ = sum(q_m_nj.*repmat(shiftdim(reflCoeffs,-2),[nFilters,nChannels,1]),3);
        minimize sum(sum_square(measVals - E*sum(measVals - coeffQ) - coeffQ)) +  alpha*sum(sum_square(R*basisRefl*reflCoeffs))
        subject to
            sum(E) == 1
            E >= 0          % This constraint is not in the paper, but we add it nonetheless.
    cvx_end
    
    hist.objValsChr(i) = cvx_optval;
    
    if abs(hist.objValsChr(i) - hist.objValsRefl(i)) < inputs.eps, break; end;

end
hist.objValsRefl = hist.objValsRefl(1:i);
hist.objValsChr = hist.objValsChr(1:i);

reflEst = basisRefl*reflCoeffs;


end

