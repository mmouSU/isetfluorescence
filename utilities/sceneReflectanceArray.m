function [scene, rcSize] = sceneReflectanceArray(reflectance,pSize,wave)

% Create a scene from an array of surface reflectance spectra. The
% reflectance input is a 3D array (h x w x nWaves) that for every spatial
% location defines the reflectance vector. The scene can be up-sampled by
% defining pSize > 1.
%
% Copyright, VISTASOFT 2016

if ieNotDefined('pSize'),    pSize = 1; end

% Default scene
scene = sceneCreate('default',1,wave);
if ieNotDefined('wave'), wave = sceneGet(scene,'wave');
else                     scene = sceneSet(scene,'wave',wave);
end

nWave = length(wave);
defaultLuminance = 100;  % cd/m2


% Spatial arrangement
r = size(reflectance,1); c = size(reflectance,2);
rcSize = [r,c];

% Convert the scene reflectances into photons assuming an equal energy
% illuminant.
ee         = ones(nWave,1);           % Equal energy vector
e2pFactors = Energy2Quanta(wave,ee);  % Energy to photon factor

% Illuminant
illuminantPhotons = diag(e2pFactors)*ones(nWave,1);

% Convert the reflectances into photons
% Data from first file are in the left columns, second file next set of
% cols, and so forth. There may be a gray strip at the end.
% Scale reflectances by incorporating energy to photon scale factr
radiance = diag(e2pFactors)*reshape(reflectance,[r*c, nWave])';    

sData = reshape(radiance',[r c nWave]);

% Store the photon data as XYZ, too. We will need these for evaluation of
% color algorithms later (sensor and illuminant correction).  These are
% stored in RGB format here (row,col,wave).
XYZ = ieXYZFromPhotons(sData,wave);

% Build up the size of the image regions - still reflectances
sData = imageIncreaseImageRGBSize(sData,pSize);

% Add data to scene, using equal energy illuminant
scene = sceneSet(scene,'photons',sData);
scene = sceneSet(scene,'illuminantPhotons',illuminantPhotons);
scene = sceneSet(scene,'illuminantComment','Equal energy');
scene = sceneSet(scene,'name','Reflectance Chart (EE)');
% vcAddObject(scene); sceneWindow;

% Adjust the illuminance to a default level in cd/m2
scene = sceneAdjustLuminance(scene,defaultLuminance);
% vcAddObject(scene); sceneWindow;

% Attach the chart parameters to the scene object so we can easily find the
% centers later
chartP.sFiles   = [];
chartP.sSamples = 1;
chartP.grayFlag = 0;
chartP.sampling = 1;
chartP.rowcol   = rcSize;
chartP.XYZ      = XYZ;
scene = sceneSet(scene,'chart parameters',chartP);

return



