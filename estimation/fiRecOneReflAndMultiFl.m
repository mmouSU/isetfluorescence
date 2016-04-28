function [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMat, predRefl, predFl, hist  ] = fiRecOneReflAndMultiFl( measVals, cameraMat, illuminant, cameraGain, basisRefl, basisEm, basisEx, alpha, beta, gamma, eta, varargin )

p = inputParser;
p.addParamValue('maxIter',200);
p.addParamValue('tauIncr',5,@isscalar);
p.addParamValue('tauDecr',5,@isscalar);
p.addParamValue('rhoInit',0.1,@isscalar);
p.addParamValue('rescaleRho',true);
p.addParamValue('mu',10,@isscalar);
p.addParamValue('epsilon',1e-6);
p.addParamValue('nFluorophores',0);
p.addParamValue('dMatRef',[]);
p.addParamValue('reflRef',[]);
p.parse(varargin{:});
inputs = p.Results;

nEmBasis = size(basisEm,2);
nExBasis = size(basisEx,2);
nReflBasis = size(basisRefl,2);
nWaves = size(illuminant,1);

%
scaleFac = max(cameraGain(:));
cameraGain = cameraGain/scaleFac;
illuminant = illuminant*scaleFac;


% Initialize variables
y1 = zeros(nWaves,1);                  % Reflectance
y2 = zeros(nWaves);                    % Fluorescence excitation-emission matrix
y3 = zeros(nEmBasis,nExBasis);         % Fluorescence excitation-emission wieghts

y1minus = y1;
y2minus = y2;
y3minus = y3;

% Dual variables
u1 = zeros(nWaves,1);
u2 = zeros(nWaves);
u3 = zeros(nEmBasis,nExBasis); 

% Estimate
sol = [];

% Pre-compute useful matrices
% [~, ~, AtA1, Atb1] = generateLSMatrices(measVals,cameraMat,illuminant,cameraGain,basisRefl,basisEx,basisEm, alpha, beta, gamma);
lowerTr = logical(tril(ones(nWaves),-1));

hist.rho = ones(inputs.maxIter+1,1);
hist.prRes = zeros(inputs.maxIter+1,1);
hist.dualRes = zeros(inputs.maxIter+1,1);
hist.pcgIter = zeros(inputs.maxIter,1);
hist.pcgAcc = zeros(inputs.maxIter,1);

hist.pixelErr = zeros(inputs.maxIter,1);

if ~isempty(inputs.dMatRef);
    hist.dMatErr = zeros(inputs.maxIter,1);
end

if ~isempty(inputs.reflRef);
    hist.reflErr = zeros(inputs.maxIter,1);
end


tElapsed = 0; t2 = tic;
for i=1:inputs.maxIter
    
    % Solve for reflectance and fluorescence
    
    % t1 = tic;
    % [~, ~, AtA2, Atb2] = generatePenaltyMatrices(y1-u1,y2-u2,y3-u3,basisRefl,basisEx,basisEm);
    
    % The iterative pcg method is 10 times faster than using function
    % handles
    % [sol, ~, hist.pcgAcc(i), hist.pcgIter(i)] = pcg(AtA1 + (hist.rho(i)/2)*AtA2,Atb1 + (hist.rho(i)/2)*Atb2,[],10000,[],[],sol);
    % tElapsed = tElapsed + toc(t1);
    
    
    
    t1 = tic;
    b = applyAtb(measVals, cameraMat, illuminant, cameraGain, basisRefl, basisEx, basisEm, y1-u1, y2-u2, y3-u3, hist.rho(i)/2);
    AtAhndl = @(x) applyAtA(x, cameraMat, illuminant, cameraGain, basisRefl, basisEx, basisEm, hist.rho(i)/2, alpha, beta, gamma);
    
    [sol, ~, hist.pcgAcc(i), hist.pcgIter(i)] = pcg(AtAhndl,b,[],10000,[],[],sol);
    tElapsed = tElapsed + toc(t1);
    
    
    reWEst = sol(1:nReflBasis);
    wEst = reshape(sol(nReflBasis+1:end),nEmBasis,nExBasis);
    
    % Solve for penalties
    % Box penalty on reflectance
    y1 = basisRefl*reWEst + u1;
    y1(y1 < 0) = 0;
    y1(y1 > 1) = 1;
    
    % Nonnegative penalty on fluorescence (i.e on the lower triangular part
    % only, we dont care about the upper part since we nuke it anyway)
    y2 = basisEm*wEst*basisEx' + u2;
    %y2 = tril(y2,-1);
    y2(lowerTr & y2 < 0) = 0;
    % y2 = sign(y2).*max(abs(y2)-(0.0002)/hist.rho(i),0);

    % If we want to strictly enforce lower triangularity:
    % (This may be too severe a limitation)
    
    % Rank penalty on fluorescence weights
    tmp = wEst+u3;
    [U,S,V] = svd(tmp);
    if ~isvector(S)
        s = diag(S);
    else
        s = S(1);
    end
    if inputs.nFluorophores <= 0
        s = sign(s).*max(abs(s)-(eta)/hist.rho(i),0);
    else
        s(inputs.nFluorophores+1:end) = 0;
    end
    % If the number of excitation end emission bases is different, then we
    % need to pad the S matrix accordingly
    y3 = U*padarray(diag(s),max([nEmBasis nExBasis] - length(s),0),0,'post')*V';
    
    % Update the scaled dual variables
    res1 = basisRefl*reWEst - y1;
    res2 = basisEm*wEst*basisEx' - y2;
    res3 = wEst - y3;
    u1 = u1 + res1;
    u2 = u2 + res2;
    u3 = u3 + res3;
    
    
    hist.prRes(i) = sqrt(norm(res1,'fro')^2 + norm(res2,'fro')^2 + norm(res3,'fro')^2);
    hist.dualRes(i) = sqrt(norm(y1-y1minus,'fro')^2 + norm(y2-y2minus,'fro')^2 + norm(y3-y3minus,'fro')^2);
    
    if (mod(i,10) == 0)
        
        tTotal = toc(t2);
        fprintf('Iter %i pr. res %f, dual res %f\n',i,hist.prRes(i),hist.dualRes(i));
        fprintf('     -> Total time: %.3f sec\n',tTotal);
        fprintf('     -> PCG time:   %.3f sec\n',tElapsed);

        t2 = tic;
        tElapsed = 0;
    end
    
    if max(hist.prRes(i),hist.dualRes(i)) < inputs.epsilon && (inputs.epsilon >=0)
        break;
    end
    
    % Re-scale the parameter rho
    if hist.prRes(i) > inputs.mu*hist.dualRes(i) && inputs.rescaleRho
        hist.rho(i+1) = hist.rho(i)*inputs.tauIncr;
    end;
    if hist.dualRes(i) > inputs.mu*hist.prRes(i) && inputs.rescaleRho
        hist.rho(i+1) = hist.rho(i)/inputs.tauDecr;
    end;
    
    u1 = u1*hist.rho(i)/hist.rho(i+1);
    u2 = u2*hist.rho(i)/hist.rho(i+1);
    u3 = u3*hist.rho(i)/hist.rho(i+1);

    
    y1minus = y1;
    y2minus = y2;
    y3minus = y3;
    
    
    %% Compute approximation error
    % The opitmization is done, let's compute the estimates
    rfCoeffs = reWEst;
    reflEst = basisRefl*rfCoeffs;
        
    dMat = max(tril(basisEm*wEst*basisEx',-1),0);

    % Compute the prediction
    predRefl = cameraGain.*(cameraMat'*diag(reflEst)*illuminant);
    predFl = cameraGain.*(cameraMat'*dMat*illuminant);
    
    
    if ~isempty(inputs.dMatRef)
        hist.dMatErr(i) = fiComputeError(dMat(:),inputs.dMatRef(:),'normalized');
    end
    
    if ~isempty(inputs.reflRef)
        hist.reflErr(i) = fiComputeError(reflEst,inputs.reflRef,'');
    end
    
    
    hist.pixelErr(i) = fiComputeError(predRefl(:) + predFl(:),measVals(:),'');
    
    
end

% The opitmization is done, let's compute the estimates
rfCoeffs = reWEst;
reflEst = basisRefl*rfCoeffs;

% Fluorescence
% W = y3;
%[U,S,V] = svd(W);

% We are returning the coefficients of the largest component
% emCoeffs = U(:,1)*sqrt(S(1,1));
%exCoeffs = V(:,1)*sqrt(S(1,1));

% Fluorescence emission and excitation shapes
% emEst = abs(basisEm*emCoeffs);
% exEst = abs(basisEx*exCoeffs);

dMat = max(tril(basisEm*wEst*basisEx',-1),0);
[U,S,V] = svd(dMat);

emEst = U(:,1)*sqrt(S(1,1))*sign(min(U(:,1)));
exEst = V(:,1)*sqrt(S(1,1))*sign(min(V(:,1)));

nF = max(exEst);
exEst = exEst/nF;
emEst = emEst*nF;

emCoeffs = basisEm'*emEst;
exCoeffs = basisEx'*exEst;

% Compute the prediction
predRefl = cameraGain.*(cameraMat'*diag(reflEst)*illuminant);

predFl = cameraGain.*(cameraMat'*dMat*illuminant);


end

function [A, b, AtA, Atb] = generateLSMatrices(measVals, camera, illuminant, cameraGain, basisRefl, basisEx, basisEm, alpha, beta, gamma)

nWaves = size(camera,1);
nFilters = size(camera,2);
nChannels = size(illuminant,2);
nReflBasis = size(basisRefl,2);

nExBasis = size(basisEx,2);
nEmBasis = size(basisEm,2);

b = measVals(:);

% Submatrix representing reflectance
rSubMat = zeros(nChannels*nFilters,nReflBasis);
for i=1:nReflBasis
    tmp = cameraGain.*(camera'*diag(basisRefl(:,i))*illuminant);
    rSubMat(:,i) = tmp(:);
end

% Submatrix representing fluorescence
cmEm = camera'*basisEm;
exLi = basisEx'*illuminant;
flSubMat = kron(exLi',cmEm);
flSubMat = diag(cameraGain(:))*flSubMat;

% Now include smoothness penalties
R = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];
% R = -[eye(nWaves - 2) zeros(nWaves-2,2)] + [zeros(nWaves-2,1) 2*eye(nWaves-2) zeros(nWaves-2,1)] -[eye(nWaves - 2) zeros(nWaves-2,2)];
rPen = R*basisRefl;

emPen = kron(basisEx,(R*basisEm));
exPen = kron(R*basisEx,basisEm);

A = [rSubMat, flSubMat;
     sqrt(alpha)*rPen, zeros(size(R,1),nExBasis*nEmBasis);
     zeros(size(emPen,1),nReflBasis) sqrt(beta)*emPen;
     zeros(size(exPen,1),nReflBasis) sqrt(gamma)*exPen];
b = [b; zeros(size(R,1) + size(emPen,1) + size(exPen,1),1)];

AtA = A'*A;
Atb = A'*b;

end


function [A, b, AtA, Atb] = generatePenaltyMatrices(reflPen,flMatPen,flWghtsPen,basisRefl,basisEx,basisEm)

A = [];
b = [];

nWaves = size(basisRefl,1);
nReflBasis = size(basisRefl,2);
nExBasis = size(basisEx,2);
nEmBasis = size(basisEm,2);

% Excitation-emission matrix
tmp = kron(basisEx,basisEm);

%{
% Reflectance
b = reflPen(:);
A = [basisRefl zeros(nWaves,nExBasis*nEmBasis)];

b = [b; flMatPen(:)];
A = [A; zeros(nWaves^2,nReflBasis) tmp];

% Excitation-emission matrix weights
b = [b; flWghtsPen(:)];
A = [A; zeros(nExBasis*nEmBasis,nReflBasis) eye(nExBasis*nEmBasis)];
%}



AtA = [basisRefl'*basisRefl zeros(nReflBasis,nExBasis*nEmBasis);
       zeros(nExBasis*nEmBasis,nReflBasis) tmp'*tmp + eye(nExBasis*nEmBasis)];
Atb = [basisRefl'*reflPen(:);
       tmp'*flMatPen(:) + flWghtsPen(:)];

end



function res = applyAtA(x, cameraMat, illuminant, cameraGain, reflBasis, exBasis, emBasis, rho, alpha, beta, gamma)

    nWaves = size(cameraMat,1);

    nReflBasis = size(reflBasis,2);
    nExBasis = size(exBasis,2);
    nEmBasis = size(emBasis,2);
    
    nabla = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];

    cameraGain2 = cameraGain.^2;
    
    wr = x(1:nReflBasis);
    W = x(nReflBasis+1:end);
    W = reshape(W,nEmBasis,nExBasis);
    
    tmp1 = reflBasis'*diag(cameraMat*(cameraGain2.*(cameraMat'*diag(reflBasis*wr)*illuminant)*illuminant'));
    tmp2 = reflBasis'*diag(cameraMat*(cameraGain2.*(cameraMat'*tril(emBasis*W*exBasis',-1)*illuminant)*illuminant'));
    
    % OLD line read: res1 = 2*tmp1 + ...
    % Empirically checked it was wrong, should go through the derivation...
    res1 = tmp1 + alpha*reflBasis'*(nabla'*nabla)*reflBasis*wr + rho*(reflBasis'*reflBasis)*wr + tmp2;
    
    tmp3 = emBasis'*tril(cameraMat*(cameraGain2.*(cameraMat'*diag(reflBasis*wr)*illuminant)*illuminant'),-1)*exBasis;
    tmp4 = emBasis'*tril(cameraMat*(cameraGain2.*(cameraMat'*tril(emBasis*W*exBasis',-1)*illuminant)*illuminant'),-1)*exBasis;
    
    %res2 = tmp3 + tmp4 + beta*emBasis'*(nabla'*nabla)*emBasis*W*(exBasis'*exBasis) + ...
    %                   + gamma*(emBasis'*emBasis)*W*exBasis'*(nabla'*nabla)*exBasis + ...
    %                   + rho*(emBasis'*emBasis)*W*(exBasis'*exBasis) + ...
    %                   + rho*W;
                   
    res2 = tmp3 + tmp4 + beta*emBasis'*tril((nabla'*nabla)*tril(emBasis*W*exBasis',-1),-1)*exBasis + ...
                       + gamma*emBasis'*tril(tril(emBasis*W*exBasis',-1)*(nabla'*nabla),-1)*exBasis + ...
                       + rho*(emBasis'*emBasis)*W*(exBasis'*exBasis) + ...
                       + rho*W;
                   
    res = [res1(:); res2(:)];
    
end

function res = applyAtb(measVals, cameraMat, illuminant, cameraGain, reflBasis, exBasis, emBasis, y1, y2, y3, rho)

    res1 = reflBasis'*diag(cameraMat*(cameraGain.*measVals)*illuminant') + rho*reflBasis'*y1;
    
    res2 = emBasis'*tril(cameraMat*(cameraGain.*measVals)*illuminant',-1)*exBasis + ...
           rho*emBasis'*y2*exBasis + rho*y3;
       
    res = [res1(:); res2(:)];

end
