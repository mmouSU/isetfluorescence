function [ refl, fl ] = fiComputeReflFlContrib( camera, illuminant, cameraGain, reflectance, dm )

nSamples = size(reflectance,2);
nFilters = size(camera,2);
nChannels = size(illuminant,2);

refl = zeros(nFilters, nChannels, nSamples);
fl = zeros(nFilters, nChannels, nSamples);

for s=1:nSamples

    refl(:,:,s) = cameraGain.*(camera'*diag(reflectance(:,s))*illuminant);
    fl(:,:,s) = cameraGain.*(camera'*dm{s}*illuminant);

end


end

