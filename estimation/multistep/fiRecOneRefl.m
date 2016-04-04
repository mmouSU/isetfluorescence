function [ reflEst, reflCoeffs, E, hist ] = fiRecOneRefl( measVals, cameraMat, cameraGain, illuminant, basisRefl, alpha, varargin )

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
    
    cvx_begin quiet
        variable reflCoeffs(nReflBasis,1)
        coeffQ = sum(q_m_nj.*repmat(shiftdim(reflCoeffs,-2),[nFilters,nChannels,1]),3);
        minimize sum(sum_square(measVals - E*sum(measVals - coeffQ) - coeffQ)) +  alpha*sum(sum_square(R*basisRefl*reflCoeffs))
        subject to
            basisRefl*reflCoeffs >= 0
    cvx_end
    
    hist.objValsRefl(i) = cvx_optval;
    
    cvx_begin quiet
        variable E(nFilters,1)
        coeffQ = sum(q_m_nj.*repmat(shiftdim(reflCoeffs,-2),[nFilters,nChannels,1]),3);
        minimize sum(sum_square(measVals - E*sum(measVals - coeffQ) - coeffQ)) +  alpha*sum(sum_square(R*basisRefl*reflCoeffs))
        subject to
            sum(E) == 1
            E >= 0          % This constraint is not discussed in the paper.
    cvx_end
    
    hist.objValsChr(i) = cvx_optval;
    
    if abs(hist.objValsChr(i) - hist.objValsRefl(i)) < inputs.eps, break; end;

end
hist.objValsRefl = hist.objValsRefl(1:i);
hist.objValsChr = hist.objValsChr(1:i);

reflEst = basisRefl*reflCoeffs;


end

