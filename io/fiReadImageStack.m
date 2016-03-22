function [ RAW, cameraGain, scaledRAW  ] = fiReadImageStack( fName )

data = load(fName);

h = size(data.RAW,1);
w = size(data.RAW,2);

RAW = im2double(data.RAW)/255;

cameraGain = data.shutter.*(10.^(data.gain/20));

scaledRAW = RAW./repmat(shiftdim(cameraGain,-2),[h w 1 1]);


end

