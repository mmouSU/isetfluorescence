% s_TeethFluorescenceData.m

% Before we put in a filter to block the reflected light, the measurements
% of the radiance from tongue and teeth under the excitation light had two
% components.  One component was reflectance (the amount of light reflected from the
% tongue (or teeth) under the excitation light) and one component fluorescence (the amount of light emitted from
% the tongue (or teeth) under the excitation light. 
% In order to get an estimate of fluorescence, we need to subtract out the
% component due to reflectance.
% The way we do this is to first measure tongue reflectance and then
% multiply tongue reflectance by the excitation light. This will be the
% amount of light we expect to be reflected when the tongue is illuminated
% with the excitation light.
% Finally, we measure the radiance of the tongue under the excitation light
% and subtract out the expected reflectance component.  

% This is what we do in this script. 

% Note that there will be some amount of error due to the fact that the
% tongue might not be in the same position for both tungsten and the
% excitation light and also not in the same position as the white surface
% we used to measure the spectral energy of the two lights. 

cd /users/joyce/Github/isetcam/;
addpath(genpath(pwd));
cd /users/joyce/Github/oraleye;
addpath(genpath(pwd));
cd /users/joyce/Github/isetfluorescence;
addpath(genpath(pwd));

wave = (350:5:700);
OralEyeLight = ieReadSpectra('OralEye_385.mat',wave);

%% Subject 001 

fname = fullfile(fiToolboxRootPath,'data','teeth','Subject001','TungstenLight');
TungstenLight= ieReadSpectra(fname,wave);

fname = fullfile(fiToolboxRootPath,'data','teeth','Subject001',' BlueFlashlight');
BlueFlashlight = ieReadSpectra(fname,wave);

% Compare the lights
% ieNewGraphWin; plot(wave,OralEyeLight); hold on;
% plot(wave, BlueFlashlight);

% Teeth
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject001',' TeethRadianceUnderTungsten');
TeethRadianceUnderTungsten = ieReadSpectra(fname,wave);

TeethReflectance = TeethRadianceUnderTungsten ./ TungstenLight;
% ieNewGraphWin; plot(wave,TeethReflectance);

%   calculate ReflectedExcitationLight
%       the amount of light that would be reflected from excitation
%       light by multiplying the excitation light with the calculated reflectance
ReflectedExcitation = BlueFlashlight .* TeethReflectance;

%   read in TeethRadiance under excitation light
fname= fullfile(fiToolboxRootPath,'data','teeth','Subject001',' TeethRadianceUnderBlueFlashlight');
TeethRadianceUnderBlueFlashlight = ieReadSpectra(fname,wave);

%   calculate TeethFluorescence
%       the amount of fluorescence by subtracting the
%       ReflectedExcitationLight from TeethRadiance 
TeethFluorescence = TeethRadianceUnderBlueFlashlight - ReflectedExcitation;
ieNewGraphWin; plot(wave,TeethFluorescence/max(TeethFluorescence),'r','linewidth',2); axis([380 700 0.0 1]); hold on;


%% Subject 007 
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject007','TungstenLight');
TungstenLight= ieReadSpectra(fname,wave);

fname = fullfile(fiToolboxRootPath,'data','teeth','Subject007',' BlueFlashlight');
BlueFlashlight = ieReadSpectra(fname,wave);

% Compare the lights
% ieNewGraphWin; plot(wave,OralEyeLight); hold on;
% plot(wave, BlueFlashlight);

% Teeth
%   read in TeethReflectance under tungsten light and calculate reflectance
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject007',' TeethRadianceUnderTungsten');
TeethRadianceUnderTungsten = ieReadSpectra(fname,wave);

TeethReflectance = TeethRadianceUnderTungsten ./ TungstenLight;
% ieNewGraphWin; plot(wave,TeethReflectance);

%   calculate ReflectedExcitationLight
%       the amount of light that would be reflected from excitation
%       light by multiplying the excitation light with the calculated reflectance
ReflectedExcitation = BlueFlashlight .* TeethReflectance;

%   read in TeethRadiance under excitation light
fname= fullfile(fiToolboxRootPath,'data','teeth','Subject007',' TeethRadianceUnderBlueFlashlight');
TeethRadianceUnderBlueFlashlight = ieReadSpectra(fname,wave);

%   calculate TeethFluorescence
%       the amount of fluorescence by subtracting the
%       ReflectedExcitationLight from TeethRadiance 
TeethFluorescence = TeethRadianceUnderBlueFlashlight - ReflectedExcitation;
plot(wave,TeethFluorescence/max(TeethFluorescence),'g','linewidth',2); axis([380 700 0.0 1]); hold on;


%% Subject 008 
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject008','TungstenLight');
TungstenLight= ieReadSpectra(fname,wave);

fname = fullfile(fiToolboxRootPath,'data','teeth','Subject008',' BlueFlashlight');
BlueFlashlight = ieReadSpectra(fname,wave);

% Compare the lights
% ieNewGraphWin; plot(wave,OralEyeLight); hold on;
% plot(wave, BlueFlashlight);

% Teeth
%   read in TeethReflectance under tungsten light and calculate reflectance
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject008',' TeethRadianceUnderTungsten');
TeethRadianceUnderTungsten = ieReadSpectra(fname,wave);

TeethReflectance = TeethRadianceUnderTungsten ./ TungstenLight;
% ieNewGraphWin; plot(wave,TeethReflectance);

%   calculate ReflectedExcitationLight
%       the amount of light that would be reflected from excitation
%       light by multiplying the excitation light with the calculated reflectance
ReflectedExcitation = BlueFlashlight .* TeethReflectance;

%   read in TeethRadiance under excitation light

fname= fullfile(fiToolboxRootPath,'data','teeth','Subject008',' TeethRadianceUnderBlueFlashlight');
TeethRadianceUnderBlueFlashlight = ieReadSpectra(fname,wave);

%   calculate TeethFluorescence
%       the amount of fluorescence by subtracting the
%       ReflectedExcitationLight from TeethRadiance 
TeethFluorescence = TeethRadianceUnderBlueFlashlight - ReflectedExcitation;
plot(wave,TeethFluorescence/max(TeethFluorescence),'b','linewidth',2); axis([380 700 0.0 1]); hold on;

