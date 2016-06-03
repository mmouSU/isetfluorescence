% Evaluate the single fluorophore algorithm convergence using simulated
% data and a model of the real acquisition system.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Save results to file if saveFName ~= []
% saveFName = fullfile(fiToolboxRootPath,'results','evaluation','conv_Fl.mat');
saveFName = [];

inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

deltaL = wave(2) - wave(1);
nWaves = length(wave);

alpha = 0.01;
beta = 0.01;

maxIter = 1000;                      % Set the maximal number of algorithm iterations

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

[ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, reflValsEst, flValsEst, hist  ] = ...
fiRecReflAndFl( measVals, camera, cameraGain*deltaL, cameraOffset, illuminant, reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter',maxIter,...
'reflRef',reflRef,'exRef',exRef,'emRef',emRef,'eps',0,'pixelRef',true);

if ~isempty(saveFName)
    save(saveFName,'reflEst','emEst','exEst','reflValsEst','flValsEst','hist',...
        'maxIter','alpha','beta','exRef','emRef','reflRef','inFName');
end
       
      