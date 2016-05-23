function [ emEst ] = fiRecOneEm( flChromaticity, cameraMat, DB )

% [ emEst ] = fiRecOneEm( flChromaticity, cameraMat, DB )
%
% Search a databse of fluorescence emission spectra and find those that
% have the most similar chromaticities with the input. This function
% implements one of the steps of the estimation algorithm by Fu et al. CVPR
% 2014, section 3.4
%
% Inputs:
%    flChromaticity - a (c x s) array of fluorescence emission
%      chromaticities of s different surfaces.
%    cameraMat - a (w x c) matrix containing the spectral responsivity
%      functions of the c camera channels sampled at w wavebands.
%    DB - a (w x t) array containing a collection of t emission spectra. The 
%      emission spectrum estimate will be selected from this data set.
%
% Outputs:
%    emEst - a (w x s) array of fluorescence emission spectra for the s
%      surfaces.
%
% Copyright, Henryk Blasinski 2014.

dbSize = size(DB,2);

estChr = cameraMat'*DB;
estChr = estChr*diag(1./sum(estChr));

% Compute the distance in the chromaticity space.
dist = sum((estChr - repmat(flChromaticity,[1 dbSize])).^2);

% Pick the smallest distance.
[~, pos] = min(dist);
emEst = DB(:,pos);

end

