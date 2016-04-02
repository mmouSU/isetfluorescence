function [ RAW, cameraGain, scaledRAW, shutter, gain  ] = fiReadImageStack( fName )

data = load(fName);

h = size(data.RAW,1);
w = size(data.RAW,2);

RAW = im2double(data.RAW)/255;
shutter = data.shutter/1000; % In seconds
gain = data.gain;

cameraGain = shutter.*(10.^(gain/20));

scaledRAW = RAW./repmat(shiftdim(cameraGain,-2),[h w 1 1]);


end

