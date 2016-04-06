close all;
clear variables;
clc;

nWaves = 10;
nChannels = 10;
nFilters = 10;
alpha = 0.5;
sigma = 0.0000;

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
        
figure;
hold on; grid on; box on;
plot(R(:),Rcvx(:),'.');
xlabel('True');
ylabel('cvx estimate');
title('Reflectance');

figure;
hold on; grid on; box on;
plot(F(:),Fcvx(:),'.');      
xlabel('True');
ylabel('cvx estimate');
title('Fluorescence');

%% Nuclear norm minimization

[ reflectanceEst, emissionEst, excitationEst, donaldsonMat, predRefl, predFl, hist  ] = ...
    recoverOneReflAndFlNucNorm( measVals, cameraMat, cameraGain, illuminant, alpha, sigma,...
    'maxIter',10000,...
    'updateRho',false,...
    'rhoInit',100);

figure;
hold on; grid on; box on;
plot(Rcvx(:),reflectanceEst(:),'.');
xlabel('cvx estimate');
ylabel('admm estimate');
title('Reflectance');

figure;
hold on; grid on; box on;
plot(Fcvx(:),donaldsonMat(:),'.');
xlabel('cvx estimate');
ylabel('admm estimate');
title('Fluorescence');

figure;
hold on; grid on; box on;
plot([hist.prRes hist.dualRes]);
