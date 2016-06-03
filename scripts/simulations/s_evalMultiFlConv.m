% Evaluate the multi-fluorophore algorithm convergence using simulated
% data and a model of the real acquisition system.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Save the results to file if the saveFName is not empty
% saveFName = fullfile(fiToolboxRootPath,'results','evaluation','conv_multiFl.mat');
saveFName = [];

inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

% Tuning parameters
alpha = 0.01;
beta = 0.01;
eta = 0.01;

maxIter = 1000;

% Create basis function sets
nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = fiCreateBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = fiCreateBasisSet('excitation','wave',wave','n',nExBasis);
[emBasis, emScore] = fiCreateBasisSet('emission','wave',wave','n',nEmBasis);


% Load the light spectra (in photons)
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
illuminant = Energy2Quanta(wave,illuminant);
nChannels = size(illuminant,2);

% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);
       
%% Load simulation data 

nSamples = size(measVals,3);
cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

[ reflEst, reflCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMatEst, reflValsEst, flValsEst, hist  ] = ...
    fiRecReflAndMultiFl( measVals, camera, illuminant, cameraGain*deltaL,...
                         cameraOffset, reflBasis, emBasis, exBasis, alpha, beta, beta, eta, 'maxIter',maxIter,...
                         'dMatRef',dMatRef,'reflRef',reflRef,'pixelRef',true,'rescaleRho',false);

%% Save results
if ~isempty(saveFName)
    save(saveFName,'reflEst','emEst','exEst','dMatEst','reflValsEst','flValsEst','hist',...
        'maxIter','alpha','beta','eta','dMatRef','reflRef','inFName');
end


