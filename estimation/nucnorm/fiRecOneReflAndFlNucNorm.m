function [ R, emissionEst, excitationEst, F, predRefl, predFl, hist  ] = fiRecOneReflAndFlNucNorm( measVals, cameraMat, cameraGain, illuminant, alpha, sigma, varargin )

% [ R, emissionEst, excitationEst, F, predRefl, predFl, hist  ] = fiRecOneReflAndFlNucNorm( measVals, cameraMat, cameraGain, illuminant, alpha, sigma, ... )
%
% Estimate the reflectance and fluorescence of a single surface using an
% adaptation of the algorithm from Suo et al. 'Bispectral coding:
% compressive and high0quality acquisition of fluorescence and
% reflectance,' Optics Express 2014. In this implementation we solve for
% spectral properties (F and R) of a single pixel, i.e. we solve the
% optimization problem (5) using our own ADMM solver. Our attemnt was to
% follow this paper to the letter. Specifically we do not impose any
% additional constratins that could follow from the physics of the problem.
% For example the matrix R is not restricted to being diagonal nor
% non-negative.
%
% Inputs:
%   measVals - a (f x c) array of pixel measurements
%   cameraMat - a (w x f) array of spectral responsivity
%      functions of the f camera filters.
%   cameraGain - a (f x c) array of camera gains applied to
%      every channel-filter combination.
%   illuminant - a (w x c) array of the illuminant spectral
%      power distributions of the c illuminants.
%   alpha - a scalar controling the l1 norm applied to the reflectance
%      matrix R.
%   sigma - a scalar controlling the admissible error in the predicted
%      pixel values.
%
% Inputs (optional):
%   'maxIter' - the maximal number of ADMM interations (default = 500).
%   'rescaleRho' - implement a heuristic algorithm improving ADMM
%      convergence from Boyd 2011 (default = true).
%   'tauIncr' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 2).
%   'tauDecr' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 2).
%   'rhoInit' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 1).
%   'mu' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 10).
%   'verbose' - display residual values every 10th iteration 
%      (default = true).
%   'epsilon' - the convergence threshold for the ADMM algorithm 
%      (default = 1e-6).
%
% Outputs:
%   R - an (w x w) square array representing the surface spectral
%      reflectance. Note the algorithm does not impose any constraints on
%      this array, it does not have to be diagonal or non-negative.
%   emissionEst - an (w x 1) emission spectrum estimation.
%   excitationEst - an (w x 1) excitation spectrum estimation.
%   F - a (w x w) Donaldson matrix estimate.
%   predRefl - a (f x c) array of predicted pixel
%      intensities due to reflected light component.
%   predFl - a (f x c) array of predicted pixel
%      intensities due to fluoresced light component.
%   hist - a structure containing primal and dual residuals of the ADMM
%      solver, as well as number of iterations and errors of the PCG
%      algorithm used to solve least-squares minimization at every ADMM
%      iteration.
% 
% Copyright, Henryk Blasinski 2016.


p = inputParser;
p.addParamValue('maxIter',500);
p.addParamValue('tauIncr',2,@isscalar);
p.addParamValue('tauDecr',2,@isscalar);
p.addParamValue('rhoInit',1,@isscalar);
p.addParamValue('rescaleRho',true);
p.addParamValue('mu',10,@isscalar);
p.addParamValue('verbose',true);
p.addParamValue('epsilon',1e-6);
p.parse(varargin{:});
inputs = p.Results;


% Re-scale the numbers to improve numerical stability
nF = max(illuminant(:));
illuminant = illuminant/nF;
cameraGain = cameraGain*nF;


nFilters = size(cameraMat,2);
nChannels = size(illuminant,2);
nWaves = size(illuminant,1);

Y1 = zeros(nWaves);
Y2 = zeros(nWaves);
Y3 = zeros(nFilters,nChannels);

Y1minus = zeros(nWaves);
Y2minus = zeros(nWaves);
Y3minus = zeros(nFilters,nChannels);

U1 = zeros(nWaves);
U2 = zeros(nWaves);
U3 = zeros(nFilters,nChannels);
U4 = zeros(nFilters,nChannels);

sol = [];

hist.rho = inputs.rhoInit*ones(inputs.maxIter+1,1);
hist.prRes = zeros(inputs.maxIter,1);
hist.dualRes = zeros(inputs.maxIter,1);
hist.pcgIter = zeros(inputs.maxIter,1);
hist.pcgAcc = zeros(inputs.maxIter,1);

t2 = 0; t3 = tic;
for i=1:inputs.maxIter
    
    % Primal variable (F,R,N) optimization
    t1 = tic;
    b = applyAtb(cameraMat, illuminant, cameraGain, Y1 - U1, Y2 - U2,Y3 - U3,measVals - U4);
    Ahndl = @(x) applyAtA(x, cameraMat, illuminant, cameraGain, nWaves, nFilters, nChannels);
    
    [sol, ~, hist.pcgAcc(i), hist.pcgIter(i)] = pcg(Ahndl,b,[],1000,[],[],sol);
    t2 = t2 + toc(t1);
    
    F = reshape(sol(1:nWaves^2),nWaves,nWaves);
    R = reshape(sol(nWaves^2+1:2*nWaves^2),nWaves,nWaves);
    N = reshape(sol(2*nWaves^2+1:end),nFilters,nChannels);
    
    if isnan(sum(F+U1))
        fprintf('Error');
    end
    
    % Y1 update: singular value soft thresholding 
    [Utmp, Stmp, Vtmp] = svd(F + U1);
    Y1 = Utmp*softThresholding(Stmp,1/hist.rho(i))*(Vtmp');
    
    % Y2 update: soft thresholding
    Y2 = softThresholding(R + U2,alpha/hist.rho(i));
    
    % Y3 update: projection
    Y3 = N + U3;
    Y3(Y3 < -3*sigma) = -3*sigma;
    Y3(Y3 > 3*sigma) = 3*sigma;
    
    % Variable update
    U1 = U1 + F - Y1;
    U2 = U2 + R - Y2;
    U3 = U3 + N - Y3;
    U4 = U4 + cameraGain.*(cameraMat'*(R+F)*illuminant) + N - measVals;
    
    hist.prRes(i) = sqrt(norm(F - Y1,'fro')^2 + norm(R - Y2,'fro')^2 + norm(N - Y3,'fro')^2 + norm(measVals - cameraGain.*(cameraMat'*(R+F)*illuminant) - N,'fro')^2);
    hist.dualRes(i) = sqrt(norm(Y1 - Y1minus,'fro')^2 + norm(Y2 - Y2minus,'fro')^2 + norm(Y3 - Y3minus,'fro')^2);
    
    if (inputs.verbose == true) && (mod(i,10) == 0)
        fprintf('Iter %i pr. res %f, dual res %f\n',i,hist.prRes(i),hist.dualRes(i));
        fprintf('    -> PCG: iter %i, acc %f, time %f\n',hist.pcgIter(i),hist.pcgAcc(i),t2);
        fprintf('    -> Total time: %f\n',toc(t3));
        t2 = 0;
        t3 = tic;
    end
    
    if max(hist.prRes(i),hist.dualRes(i)) <= inputs.epsilon,
        break; 
    end;
    
    % Re-scale the parameter rho
    if hist.prRes(i) > inputs.mu*hist.dualRes(i) && inputs.rescaleRho
        hist.rho(i+1) = hist.rho(i)*inputs.tauIncr;
    end;
    if hist.dualRes(i) > inputs.mu*hist.prRes(i) && inputs.rescaleRho
        hist.rho(i+1) = hist.rho(i)/inputs.tauDecr;
    end;
    
    U1 = U1*hist.rho(i)/hist.rho(i+1);
    U2 = U2*hist.rho(i)/hist.rho(i+1);
    U3 = U3*hist.rho(i)/hist.rho(i+1);
    U4 = U4*hist.rho(i)/hist.rho(i+1);
    
    
    Y1minus = Y1;
    Y2minus = Y2;
    Y3minus = Y3;
    
end

hist.rho = hist.rho(1:i);
hist.prRes = hist.prRes(1:i);
hist.dualRes = hist.dualRes(1:i);
hist.pcgIter = hist.pcgIter(1:i);
hist.pcgAcc = hist.pcgAcc(1:i);


[U, ~, V] = svd(F);
emissionEst = U(:,1);
excitationEst = V(:,1);

predRefl = cameraGain.*(cameraMat'*R*illuminant);
predFl = cameraGain.*(cameraMat'*F*illuminant);

end


%% Helper functions

function res = softThresholding(in,thr)

    res = sign(in).*max(abs(in) - thr,0);

end

% Functions used by the PCG algorithm to iteratively solve a least-squares problem
% without formulating large matrices.

function res = applyAtA(x,cameraMat,illuminant,cameraGain,nWaves,nFilters,nChannels)


    F = reshape(x(1:nWaves^2),nWaves,nWaves);
    R = reshape(x(nWaves^2+1:2*nWaves^2),nWaves,nWaves);
    N = reshape(x(2*nWaves^2+1:end),nFilters,nChannels);
    
    tmp = (cameraMat*((cameraGain.^2).*(cameraMat'*(R+F)*illuminant))*illuminant') + (cameraMat*(cameraGain.*N)*illuminant');
    
    res1 = F + tmp;
    res2 = R + tmp;
    res3 = cameraGain.*(cameraMat'*(R+F)*illuminant) + 2*N;
    
    res = [res1(:); res2(:); res3(:)];

end

function res = applyAtb(cameraMat,illuminant,cameraGain,v1,v2,v3,v4)
    
    tmp = (cameraMat*(cameraGain.*v4)*illuminant');

    res1 = v1 + tmp;
    res2 = v2 + tmp;
    res3 = v3 + v4;
    
    res = [res1(:); res2(:); res3(:)];

end

