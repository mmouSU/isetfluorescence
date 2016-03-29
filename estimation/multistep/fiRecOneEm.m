function [ emEst ] = fiRecOneEm( flChromaticity, cameraMat, DB )

dbSize = size(DB,2);

estChr = cameraMat'*DB;
estChr = estChr*diag(1./sum(estChr));

dist = sum((estChr - repmat(flChromaticity,[1 dbSize])).^2);

[~, pos] = min(dist);
emEst = DB(:,pos);

end

