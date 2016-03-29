function [exEst, exCoeffs] = fiRecOneEx(measVals, reflCoeffs, chrCoeffs, cameraMat, cameraGain, illuminant, basisRefl, basisEx, gamma)

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

