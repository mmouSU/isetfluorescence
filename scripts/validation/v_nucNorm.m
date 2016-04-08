% This script checks the correctness of the ADMM implementation of the
% nuclear norm minimization algorithm described in Suo et al.
% We start by generating synthetic and unrealistic data, and then solve for
% the unknown quantities using general purpose solver such as CVX and our 
% ADMM implementation. If the implementation is correct, the results 
% between ADMM and CVX agree.
%
% Copyright, Henryk Blasinski 2016.

close all;
clear all;
clc;

% Define problem dimensionality and parameters
nWaves = 150;
nChannels = 14;
nFilters = 8;
alpha = 0.5;
sigma = 0.0000;

% Generate random data
cameraMat = eye(nWaves,nFilters);
cameraGain = rand(nFilters,nChannels);
cameraOffset = zeros(nFilters,nChannels);
illuminant = eye(nWaves,nChannels);

R = diag(rand(nWaves,1));
F = rand(nWaves,1)*rand(1,nWaves);

measVals = cameraGain.*(cameraMat'*(R+F)*illuminant);
%% CVX

cvx_begin
    variables Fcvx(nWaves,nWaves) Rcvx(nWaves,nWaves) N(nFilters,nChannels) Y1(nWaves,nWaves) Y2(nWaves,nWaves) Y3(nFilters,nChannels)
    minimize norm_nuc(Y1) + alpha*sum(norms(Y2,1,1))
    subject to
        Fcvx - Y1 == 0
        Rcvx - Y2 == 0
        N - Y3 == 0
        measVals == cameraGain.*(cameraMat'*(Rcvx + Fcvx)*illuminant) + N + cameraOffset
        -3*sigma <= N-cameraOffset <= 3*sigma
cvx_end
        

%% ADMM

[ Radmm, ~, ~, Fadmm, ~, ~, hist  ] = ...
    fiRecOneReflAndFlNucNorm( measVals, cameraMat, cameraGain, illuminant, alpha, sigma,...
    'maxIter',3000,...
    'rescaleRho',false,...
    'rhoInit',100);

%% Compare results. 
% When the results match, the figures should display identity lines.

figure;
hold on; grid on; box on;
mVal = max([Rcvx(:); Radmm(:)]);
plot(Rcvx(:),Radmm(:),'.');
plot(linspace(0,mVal,10),linspace(0,mVal,10),'r');
xlim([0 mVal]);
ylim([0 mVal]);
xlabel('cvx estimate');
ylabel('ADMM estimate');
title('Reflectance');

figure;
hold on; grid on; box on;
mVal = max([Fcvx(:); Fadmm(:)]);
plot(Fcvx(:),Fadmm(:),'.');
plot(linspace(0,mVal,10),linspace(0,mVal,10),'r');
xlim([0 mVal]);
ylim([0 mVal]);
xlabel('cvx estimate');
ylabel('ADMM estimate');
title('Donaldson matrix');

figure;
hold all; grid on; box on;
plot([hist.prRes hist.dualRes]);
set(gca,'yscale','log');
title('Convergence');

