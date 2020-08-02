% s_TeethFluorescenceData.m
%
% Before we put in a filter to block the reflected light, the measurements
% of the radiance from tongue and teeth under the excitation light had two
% components.  One component was reflectance (the amount of light reflected
% from the tongue (or teeth) under the excitation light) and one component
% fluorescence (the amount of light emitted from the tongue (or teeth)
% under the excitation light.
%
% In order to get an estimate of fluorescence, we need to subtract out the
% component due to reflectance.
%
% The way we do this is to first measure tongue reflectance and then
% multiply tongue reflectance by the excitation light. This will be the
% amount of light we expect to be reflected when the tongue is illuminated
% with the excitation light. Finally, we measure the radiance of the tongue
% under the excitation light and subtract out the expected reflectance
% component.
%
% This is what we do in this script. 
%
% Note that there will be some amount of error due to the fact that the
% tongue might not be in the same position for both tungsten and the
% excitation light and also not in the same position as the white surface
% we used to measure the spectral energy of the two lights.
%
% Joyce Farrell swears this is true

%%
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
TeethFluorescence1 = TeethRadianceUnderBlueFlashlight - ReflectedExcitation;
ieNewGraphWin; plot(wave,TeethFluorescence1,'r','linewidth',2); hold on;

%%
TeethFluorescence1 = TeethFluorescence1/max(TeethFluorescence1);
TeethFluorescence1 = ieClip(TeethFluorescence1,0,1);
TeethFluorescence1(isnan(TeethFluorescence1)) = 0;

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
TeethFluorescence2 = TeethRadianceUnderBlueFlashlight - ReflectedExcitation;
plot(wave,TeethFluorescence2/max(TeethFluorescence2),'g','linewidth',2); axis([380 700 0.0 1]); hold on;

%%
TeethFluorescence2 = TeethFluorescence2/max(TeethFluorescence2);
TeethFluorescence2 = ieClip(TeethFluorescence2,0,1);
TeethFluorescence2(isnan(TeethFluorescence2)) = 0;

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
TeethFluorescence3 = TeethRadianceUnderBlueFlashlight - ReflectedExcitation;
plot(wave,TeethFluorescence3/max(TeethFluorescence3),'b','linewidth',2); axis([380 700 0.0 1]); hold on;

%%
TeethFluorescence3 = TeethFluorescence3/max(TeethFluorescence3);
TeethFluorescence3 = ieClip(TeethFluorescence3,0,1);
TeethFluorescence3(isnan(TeethFluorescence3)) = 0;

%% Average fluorescence
%
% We are going to use this fluorescence in all the columns of the EEM from
% 380 to 400
%

excitation = zeros(size(wave));
excitation(wave < 430) = 1;
emission = (TeethFluorescence1 + TeethFluorescence2 + TeethFluorescence3)/3;

f = fluorophoreCreate();
f = fluorophoreSet(f,'wave',wave);
f = fluorophoreSet(f,'emission photons',Energy2Quanta(wave,emission(:)));
f = fluorophoreSet(f,'excitation photons',Energy2Quanta(wave,excitation'));
f = fluorophoreSet(f,'name','Teeth measured');

% Save it

fname = fullfile(fiToolboxRootPath,'data','teeth','teeth-measured.mat');
fprintf('Saving %s\n',fname);
fluorophoreSave(fname,f,'Measured in the lab by JEF');

%% Test and Plot it

teeth = fluorophoreRead('teeth-measured');

fluorophorePlot(teeth,'donaldson matrix')
fluorophorePlot(teeth,'excitation');
fluorophorePlot(teeth,'emission photons');

%%
%% END