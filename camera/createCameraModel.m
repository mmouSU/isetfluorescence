function [sensor, optics] = createCameraModel(filterID,varargin)

p = inputParser;
p.addParamValue('wave',(400:10:1000)',@isvector);
p.addRequired('filterID',@isnumeric);

p.parse(filterID,varargin{:});
inputs = p.Results;


% Here are some of the key pixel properties of the PointGrey Flea3 camera
% with an ONSemi VITA 1300 1.3 Megapixel sensor.
wellCapacity   = 13700;                     % Electrons
conversiongain = 90*1e-6;                   % Volts/electron 
voltageSwing   = conversiongain*wellCapacity;
fillfactor     = 0.5;                       % This is a realistic, albeit made up number. 
pixelSize      = 4.8*1e-6;                  % Meters
darkvoltage    = conversiongain*4.5;        % Volts/sec
readnoise      = 0.00096;                   % Volts
rows = 1024;                                 
cols = 1280; 
dsnu =  conversiongain*30;                  % Volts (dark signal non-uniformity)
prnu = 2;                                   % Percent (ranging between 0 and 100) photodetector response non-uniformity
quantizationMethod = '8bit';
analogGain   = 1;                           % Used to adjust ISO speed
analogOffset = 0;                           % Used to account for sensor black level
noiseFlag = 2;                              % Simulate noise

% Some parameters of the lens used
focalLength = 0.07;
fNumber = 4;


% Read the calibrated sensor quantum efficiency.
fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,inputs.wave,0);


% Read filter transmissivities;
fName = fullfile(fiToolboxRootPath,'camera','filters');
transmissivities = ieReadSpectra(fName,inputs.wave,0);




sensor = sensorCreate('Monochrome');
sensor = sensorSet(sensor,'name',sprintf('Monochrome, filter %i',inputs.filterID));

% Set sensor properties
sensor = sensorSet(sensor,'wave',inputs.wave);
sensor = sensorSet(sensor,'autoExposure','on'); 
sensor = sensorSet(sensor,'noiseFlag',noiseFlag);
sensor = sensorSet(sensor,'rows',rows);
sensor = sensorSet(sensor,'cols',cols);
sensor = sensorSet(sensor,'dsnulevel',dsnu);  
sensor = sensorSet(sensor,'prnulevel',prnu); 
sensor = sensorSet(sensor,'analogGain',analogGain);     
sensor = sensorSet(sensor,'analogOffset',analogOffset);
sensor = sensorSet(sensor,'quantizationmethod',quantizationMethod);
sensor = sensorSet(sensor,'infraredfilter',ones(length(inputs.wave),1));
sensor = sensorSet(sensor,'filterspectra',transmissivities(:,inputs.filterID));

% Set pixel properties
pixel =  sensorGet(sensor,'pixel');
pixel = pixelSet(pixel,'sizesamefillfactor',[pixelSize pixelSize]);   
pixel = pixelSet(pixel,'conversiongain', conversiongain);        
pixel = pixelSet(pixel,'voltageswing',voltageSwing);                                             
pixel = pixelSet(pixel,'darkvoltage',darkvoltage) ;               
pixel = pixelSet(pixel,'readnoisevolts',readnoise);  
pixel = pixelSet(pixel,'pixelspectralqe',qe); 
sensor = sensorSet(sensor,'pixel',pixel);


sensor = pixelCenterFillPD(sensor,fillfactor);



%% Optics

oi = oiCreate;
optics = oiGet(oi,'optics');
optics = opticsSet(optics,'model','DiffractionLimited');
optics = opticsSet(optics,'off axis method','Skip');
optics = opticsSet(optics,'focallength',focalLength);
optics = opticsSet(optics,'fNumber',fNumber);

end

