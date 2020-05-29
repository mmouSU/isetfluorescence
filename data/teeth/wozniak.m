%% Teeth from Wozniak and Moore 1978
%
chdir(fullfile(fiToolboxRootPath,'data','sources','wozniak'))
wave = 200:5:800;

%% We convert excitation and emission to quanta for consistency with ISETCam

foo = load('TeethEmission');
emission = ieReadSpectra('TeethEmission',wave);

% JEF wanted to extrapolate the excitation curve linearly because the plots
% stopped abruptly.  She felt that an extrapolation was called for.
excitation = ieReadSpectra('TeethExcitation',wave);

% Set up the extrapolation for emission
lst = (wave<560);
x = wave(lst); x = x(:);
y = emission(lst); y = y(:);
emission = interp1(x,y,wave,'linear','extrap');
lst = wave > 620;
emission(lst) = 0;

% Set up the extrapolation for the excitation
lst = (wave<390);
x = wave(lst); x = x(:);
y = excitation(lst); y = y(:);

excitation = interp1(x,y,wave,'pchip','extrap');
lst = wave > 435;
excitation(lst) = 0;

%% Plot 
ieNewGraphWin; 
plot(wave,emission,'k-',wave,excitation,'r-')
grid on; xlabel('Wave (nm)'); 

%% Create the fluorophore object

f = fluorophoreCreate();
f = fluorophoreSet(f,'wave',wave);

f = fluorophoreSet(f,'emission photons',Energy2Quanta(wave,emission'));
f = fluorophoreSet(f,'excitation photons',Energy2Quanta(wave,excitation'));


f = fluorophoreSet(f,'name','Teeth measured');
%{
 fluorophorePlot(f,'excitation photons')
 fluorophorePlot(f,'emission photons')
 fluorophorePlot(f,'donaldson image');
%}

%% Save it

fname = fullfile(fiToolboxRootPath,'data','teeth','teeth.mat');
fprintf('Saving %s\n',fname);
fluorophoreSave(fname,f,'digitized by JEF from Wozniak and Moore 1978 in the Journal of Dental Research');

%% END