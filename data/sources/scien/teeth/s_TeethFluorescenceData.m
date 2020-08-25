% s_TeethFluorescenceData.m
%
%  Purpose: create an EEM for tooth fluorescence
%           This script creates an EEM based on teeth radiance measurements
%           made by Joyce Farrell on 08302018
%           The EEM is based on excitation light only, but we hope to
%           improve the EEM by incorporating measurements that were made
%           and published by others using different excitation wavelengths
%           see, for example, Sundstrom, Folke, K. Fredriksson, S. Montan, Ulrika Hafstrom-Bjorkman, and Jane Strom. 1985. 
%           “Laser-Induced Fluorescence from Sound and Carious Tooth Substance: Spectroscopic Studies.” Swedish Dental Journal 9 (2): 71–80.
%
%  We illuminated a tooth in an individual with a tungsten light to measure reflectance
%  and we illuminated the same tooth (same approximate location) with a blue flashlight to measure fluorescence
%   The measurements were made with the PR670 spectroradiometer.
%   Tooth reflectance
%       The tooth reflectance was estimated by dividing the spectral radiance reflected
%       from the tooth when it was illuminated with tungsten light with the
%       spectral radiance of the tungsten light
%   Tooth fluorescence
%       The radiance measurements of the tooth under the bluelight has two components
%       One component is the light reflected from the tooth when it is illuminated with
%       the bluelight and the other component is the light emitted from the
%       tooth when it is illuminated with the bluelight
%       To estimate fluorescence, we need to subtract out the
%       component due to reflectance.
%       We estimate the bluelight reflected from the tooth by multiplying the
%       spectral energy in the bluelight with the tooth reflectance function
%       Then, we subtract the estimated reflected light from the total spectral
%       radiance measured when the tooth was illuminated with the blue light
%
%   The spectroradiometer measurement is centered on the same tooth for the tungsten and bluelight measurements. 
%   To measure the spectral energy in the tungsten and bluelight, we placed a white calibration target in the same location as the tooth when it was measured.  
%   Some measurement noise will be due to movement of the tooth or inacurracy in the placement of the white calibration surface. 


% JEF 

%%
wave = (350:5:700);
OralEyeLight = ieReadSpectra('OralEye_385.mat',wave); % In the future, measure teeth fluorescence for this light

%% Subject 001 

fname = fullfile(fiToolboxRootPath,'data','teeth','Subject001','TungstenLight');
TungstenLight1= ieReadSpectra(fname,wave);

fname = fullfile(fiToolboxRootPath,'data','teeth','Subject001',' BlueFlashlight');
BlueFlashlight1 = ieReadSpectra(fname,wave);

% Compare the lights
% ieNewGraphWin; plot(wave,OralEyeLight); hold on;
% plot(wave, BlueFlashlight);

% Teeth
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject001',' TeethRadianceUnderTungsten');
TeethRadianceUnderTungsten = ieReadSpectra(fname,wave);

TeethReflectance1 = TeethRadianceUnderTungsten ./ TungstenLight1;
% ieNewGraphWin; plot(wave,TeethReflectance);

%   calculate ReflectedExcitationLight
%       the amount of light that would be reflected from excitation
%       light by multiplying the excitation light with the calculated reflectance
ReflectedExcitation = BlueFlashlight1 .* TeethReflectance1;

%   read in TeethRadiance under excitation light
fname= fullfile(fiToolboxRootPath,'data','teeth','Subject001',' TeethRadianceUnderBlueFlashlight');
TeethRadianceUnderBlueFlashlight = ieReadSpectra(fname,wave);

%   calculate TeethFluorescence
%       the amount of fluorescence by subtracting the
%       ReflectedExcitationLight from TeethRadiance 
TeethFluorescence1 = TeethRadianceUnderBlueFlashlight - ReflectedExcitation;
% ieNewGraphWin; plot(wave,TeethFluorescence1,'r','linewidth',2); hold on;

%%
TeethFluorescence1 = TeethFluorescence1/max(TeethFluorescence1);
TeethFluorescence1 = ieClip(TeethFluorescence1,0,1);
TeethFluorescence1(isnan(TeethFluorescence1)) = 0;
ieNewGraphWin; plot(wave,TeethFluorescence1,'r','linewidth',2); hold on;
%% Subject 007 
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject007','TungstenLight');
TungstenLight2= ieReadSpectra(fname,wave);

fname = fullfile(fiToolboxRootPath,'data','teeth','Subject007',' BlueFlashlight');
BlueFlashlight2 = ieReadSpectra(fname,wave);

% Compare the lights
% ieNewGraphWin; plot(wave,OralEyeLight); hold on;
% plot(wave, BlueFlashlight2);

% Teeth
%   read in TeethReflectance under tungsten light and calculate reflectance
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject007',' TeethRadianceUnderTungsten');
TeethRadianceUnderTungsten = ieReadSpectra(fname,wave);

TeethReflectance2 = TeethRadianceUnderTungsten ./ TungstenLight2;
% ieNewGraphWin; plot(wave,TeethReflectance);

%   calculate ReflectedExcitationLight
%       the amount of light that would be reflected from excitation
%       light by multiplying the excitation light with the calculated reflectance
ReflectedExcitation = BlueFlashlight2 .* TeethReflectance2;

%   read in TeethRadiance under excitation light
fname= fullfile(fiToolboxRootPath,'data','teeth','Subject007',' TeethRadianceUnderBlueFlashlight');
TeethRadianceUnderBlueFlashlight = ieReadSpectra(fname,wave);

%   calculate TeethFluorescence
%       the amount of fluorescence by subtracting the
%       ReflectedExcitationLight from TeethRadiance 
TeethFluorescence2 = TeethRadianceUnderBlueFlashlight - ReflectedExcitation;
% plot(wave,TeethFluorescence2/max(TeethFluorescence2),'g','linewidth',2); axis([380 700 0.0 1]); hold on;

%%
TeethFluorescence2 = TeethFluorescence2/max(TeethFluorescence2);
TeethFluorescence2 = ieClip(TeethFluorescence2,0,1);
TeethFluorescence2(isnan(TeethFluorescence2)) = 0;
plot(wave,TeethFluorescence2,'g','linewidth',2); hold on;
%% Subject 008 
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject008','TungstenLight');
TungstenLight3= ieReadSpectra(fname,wave);

fname = fullfile(fiToolboxRootPath,'data','teeth','Subject008',' BlueFlashlight');
BlueFlashlight3 = ieReadSpectra(fname,wave);

% Compare the lights
% ieNewGraphWin; plot(wave,OralEyeLight); hold on;
% plot(wave, BlueFlashlight);

% Teeth
%   read in TeethReflectance under tungsten light and calculate reflectance
fname = fullfile(fiToolboxRootPath,'data','teeth','Subject008',' TeethRadianceUnderTungsten');
TeethRadianceUnderTungsten = ieReadSpectra(fname,wave);

TeethReflectance3 = TeethRadianceUnderTungsten ./ TungstenLight3;
% ieNewGraphWin; plot(wave,TeethReflectance);

%   calculate ReflectedExcitationLight
%       the amount of light that would be reflected from excitation
%       light by multiplying the excitation light with the calculated reflectance
ReflectedExcitation = BlueFlashlight3 .* TeethReflectance3;

%   read in TeethRadiance under excitation light

fname= fullfile(fiToolboxRootPath,'data','teeth','Subject008',' TeethRadianceUnderBlueFlashlight');
TeethRadianceUnderBlueFlashlight = ieReadSpectra(fname,wave);

%   calculate TeethFluorescence
%       the amount of fluorescence by subtracting the
%       ReflectedExcitationLight from TeethRadiance 
TeethFluorescence3 = TeethRadianceUnderBlueFlashlight - ReflectedExcitation;
% plot(wave,TeethFluorescence3/max(TeethFluorescence3),'b','linewidth',2); axis([380 700 0.0 1]); hold on;

%%
TeethFluorescence3 = TeethFluorescence3/max(TeethFluorescence3);
TeethFluorescence3 = ieClip(TeethFluorescence3,0,1);
TeethFluorescence3(isnan(TeethFluorescence3)) = 0;
plot(wave,TeethFluorescence3,'b','linewidth',2); hold on;
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

%% Compare the data measured for the 3 subjects
ieNewGraphWin; plot(wave,TungstenLight1,'r','linewidth',2); hold on;
plot(wave,TungstenLight2,'g','linewidth',2);
plot(wave,TungstenLight3,'b','linewidth',2);

ieNewGraphWin; plot(wave,TeethReflectance1,'r','linewidth',2); hold on;
plot(wave,TeethReflectance2,'g','linewidth',2);
plot(wave,TeethReflectance3,'b','linewidth',2);

ieNewGraphWin; plot(wave,BlueFlashlight1,'r','linewidth',2); hold on;
plot(wave,BlueFlashlight2,'g','linewidth',2);
plot(wave,BlueFlashlight3,'b','linewidth',2);
plot(wave,OralEyeLight,'k','linewidth',2);
title('Oral Eye Light compared to the blue flashlights used for 3 subjects');

ieNewGraphWin; plot(wave,TeethFluorescence1,'r','linewidth',2); hold on;
plot(wave,TeethFluorescence2,'g','linewidth',2);
plot(wave,TeethFluorescence3,'b','linewidth',2);
legend('Subject 001 Bluelight', 'Subject 007 Bluelight', 'Subject 008 Bluelight', 'OralEye Light')



%% END