function [exEst, exCoeffs] = fiRecOneEx(measVals, reflCoeffs, chrCoeffs, cameraMat, cameraGain, illuminant, basisRefl, basisEx, gamma)

% [exEst, exCoeffs] = fiRecOneEx(measVals, reflCoeffs, chrCoeffs, cameraMat, cameraGain, illuminant, basisRefl, basisEx, gamma)
% 
% Estimate the excitation spectrum. This function implements the algorithm
% of Fu et al. CVPR 2014 paper, section 3.3.
%
% Inputs:
%    measVals - a (f x c x s) matrix containing pixel intensities of s 
%      different surfaces captured with f camera channels captured under c
%      different illuminants.
%    reflCoeffs - a (n x s) array containing the estimated reflectance
%      basis function weights for the s surfaces.
%    chrCoeffs - a (f x s) array of estimated fluorescence emission
%      chromaticities of the s surfaces.
%    cameraMat - a (w x c) matrix containing the spectral responsivity
%      functions of the c camera channels sampled at w wavebands.
%    cameraGain - a (f x c x s) matrix of linear camera model gains for each
%      filter-channel-surface combination.
%    illuminant - a (w x c) matrix containing the spectral power
%      distributions of the c illuminants. 
%    basisRefl, basisEx - (w x n) matrices of n linear basis
%      functions representing reflectance, and excitation spectra
%      respectively.
%    gamma - scalar tuning parameter controlling the smoothness of
%      fluorescence excitation estimates.
%
% Outputs:
%    exEst - a (w x s) matrix of the estimated surface excitation spectra.
%    exCoeffs - a (n x s) matrix expressing the estimated surface excitation 
%      spectra in terms of the linear basis weights.
%
% Copyright, Henryk Blasinski 2016


nExBasis = size(basisEx,2);
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

coeffQ = sum(q_m_nj.*repmat(shiftdim(reflCoeffs,-2),[nFilters,nChannels,1]),3);
fl = measVals - coeffQ;

h_m_i = illuminant'*basisEx;

R = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1) eye(nWaves-1,1)];

cvx_begin quiet
    variable exCoeffs(nExBasis,1)
    minimize sum(sum_square(fl - chrCoeffs*(h_m_i*exCoeffs)')) +  gamma*sum(sum_square(R*basisEx*exCoeffs))
        subject to
            basisEx*exCoeffs >= 0
cvx_end

exEst = basisEx*exCoeffs;



end

