function [ reflectanceEst, coeffs, predResp ] = fiRecRefl( measVals, cameraResp, cameraGain, illuminant, basisFcns, lambda )

% [ reflectanceEst, coeffs, predResp ] = fiRecRefl( measVals, cameraResp, illuminant, basisFcns, lambda )
%
% This function computes the surface reflectance estimates given camera and
% illumination parameters. It is assumed that the camera has "nFilter"
% different spectral channels and the scene is observed under "nChannel"
% illuminants. 
%
% This function estimates the reflectance by approximating the
% spectral reflectance with a small number of basis functions and computing
% a least-squares fit of these coefficients to acquited data.
%
% Inputs:
%   measVals - a nFilter x nChannel x nSamples matrix with the pixel values
%   cameraResp - a nWaves x nFilter matrix with the spectral responsivity
%   of each camera channel
%   channel-illuminant combination
%   illuminant - a nWaves x nChannel matrix with illumination intensity (in
%   photons)
%   basisFcns - a nWaves x nBasis matrix with basis function
%   representations for the spectral reflectance.
%   lambda - a scalar describing the smoothness parameter
%
% Outputs:
%   reflectanceEst - a nWaves x nSamples matrix with spectral reflectance
%   estimates
%   coeffs - a nBasis x nSamples matrix with basis function weights for
%   each reflectance estimate
%   predResp - a nFilter x nChannel x nSamples matrix with modeled pixel
%   responses.
%
% Copyright, Henryk Blasinski 2014

nBasis = size(basisFcns,2);
nChannels = size(illuminant,2);
nFilters = size(cameraResp,2);
nSamples = size(measVals,3);
nWaves = size(illuminant,1);

b = measVals(:);

% Basis vectors
basisVectors = zeros(nChannels*nFilters*nSamples,nBasis*nSamples);
for j=1:nSamples
for i=1:nBasis
    tmp = cameraGain(:,:,j).*(cameraResp'*diag(basisFcns(:,i))*illuminant);
    basisVectors((j-1)*nChannels*nFilters+1:j*nChannels*nFilters,(j-1)*nBasis + i) = tmp(:);
end
end

% Create a penalty on roughness
R = [diag(ones(nWaves-1,1)) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) diag(ones(nWaves-1,1))];
Rm = [];

for i=1:nSamples
    Rm = blkdiag(Rm,R*basisFcns);
end


% Create the image formation model matrix
A = [basisVectors; sqrt(lambda)*Rm];
b = [b; zeros((nWaves-1)*nSamples,1)];

res = A\b;
coeffs = reshape(res,nBasis,nSamples);
reflectanceEst = basisFcns*coeffs;


predResp = A(1:nFilters*nChannels*nSamples,:)*res;
predResp = reshape(predResp,[nFilters, nChannels, nSamples]);
    
end

