function [ RAW, cameraGain, scaledRAW, shutter, gain  ] = fiReadImageStack( fName )

% [ RAW, cameraGain, scaledRAW, shutter, gain  ] = fiReadImageStack( fName )
%
% Read in an image stack from a Matlab file defined in fName, saved with
% the custom acquisiton hardware. The ISET model of the hardware system is
% available in [fiToolboxRootPath]/camera
%
% Inputs:
%   fName - a path to Matlab (.mat) file containing the data
%
% Outputs:
%   RAW - a (h x w x nFilters x nChannels) array of raw sensor intensities
%      scaled to [0,1] interval. The image size is (h x w) and the scene
%      was captured under nChannels different illuminants and through
%      nFilters different bandpass filters.
%   cameraGain - a (nFilters x nChannels) array of gain factors applied by
%      the camera to every filter-channel combination. This gain factor
%      incorporates the effect of amplifier gains as well as integration
%      times.
%   scaledRAW - a (h x w x nFilters x nChannels) array of RAW data scaled by
%      corresponding cameraGain. The scaledRAW data is linear in pixel
%      responsivities, illuminant, reflectance and fluorescence.
%   shutter - a (nFilters x nChannels) array of image integration times (in
%      seconds) for every filter-channel pair.
%   gain - a (nFilters x nChannels) array of sensor amplifier gains (in dB)
%      for every filter-channel pair.
%
% Copyright, Henryk Blasinski 2016.

data = load(fName);

h = size(data.RAW,1);
w = size(data.RAW,2);

RAW = im2double(data.RAW)/255;
shutter = data.shutter/1000;    % Convert ms to seconds
gain = data.gain;               % in dB

cameraGain = shutter.*(10.^(gain/20));

scaledRAW = RAW./repmat(shiftdim(cameraGain,-2),[h w 1 1]);


end

