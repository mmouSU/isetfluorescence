%% Read in and create additional fluorophores from vectors provided by JEF
%
% These data were scanned by JEF from the Monici 2005 paper.
%
% We believe the numbers read in from the literature (always specified as
% arbitrary units) are in fact in relative energy units (not relative
% photon units).  So we convert the energy units to quanta here.
%
% Wandell

%%
chdir(fiToolboxRootPath)
chdir('local')
wave = 200:5:800;

%% We convert excitation and emission to quanta for consistency with ISETCam

% Scanned molecule names
names = {'Collagen','Elastin','Flavins','Lipopigments','NADH','Porphyrins'};

% For each type of molecule
for nn = 1:numel(names)
    fprintf('Processing %s\n',names{nn});
    
    % Emission data
    emFile = sprintf('%sEmission.mat',names{nn});
    emData = load(emFile);
    data   = emData.(names{nn});
    emWave = data(:,1);
    em     = data(:,2);
    em     = Energy2Quanta(emWave,em);
    
    % NaNs for the extrapolation region
    em = interp1(emWave,em,wave);
    em     = em/max(em(:));
    em(isnan(em)) = 0;

    % Excitation data
    % The literature thinks we have
    %    filter * energy
    % But we want a filterPhoton * Photon
    % To get that we must have filter*(unit photon/ unit energy) * energy
    % So we multiply the filter by 
    exFile = sprintf('%sExcitation.mat',names{nn});
    exData = load(exFile);
    str = sprintf('%sExcitation',names{nn});
    exWave = exData.(str)(:,1);
    ex     = exData.(str)(:,2);
    ex     = ex .* Energy2Quanta(exWave,ones(numel(exWave),1));
    % ieNewGraphWin; plot(exWave,ex);
    
    % NaNs for the extrapolation region
    ex = interp1(exWave,ex,wave);
    ex = ex/max(ex(:));
    ex(isnan(ex)) = 0;
    
    % Build the fluorophore and write it out
    f = fluorophoreCreate();
    f = fluorophoreSet(f,'wave',wave);

    f = fluorophoreSet(f,'emission photons',em);
    f = fluorophoreSet(f,'excitation photons',ex);
    
    % fluorophorePlot(f,'excitation photons')
    % fluorophorePlot(f,'emission photons')
    
    f = fluorophoreSet(f,'name',names{nn});
    % fluorophorePlot(f,'donaldson image');
    
    fname = fullfile(fiToolboxRootPath,'data','Monici',names{nn});
    
    fprintf('Saving %s\n',names{nn});
    fluorophoreSave(fname,f,'digitized by JEF from Monici 2005');
    fluorophorePlot(f,'donaldson image');

end

%%  Compare with the webfluor directly

% These are from the webfluor site
%
% The diagonal terms are the reflectance.  These are not quite zero.
%{
fname = 'FAD.txt';
fname = 'Flavin.txt';
fname = 'NADH.txt'
fname = 'NADPH.txt';
fname = 'hemoglobin.txt';
fname = 'protoporphyrin.txt';
fname = 'collagen1.txt';
fname = 'collagen4.txt';
fname = 'collagen5.txt';
fname = 'collagen6.txt';
fname = 'collagen7.txt';
fname = 'elastin.txt';
%}
% chdir('data'); chdir('webfluor');
fname = 'FAD.txt';
T = readtable(fname);
exWave = T{:,1};
emWave = 260:5:750;
exemMatrix = T{:,2:end}';
ieNewGraphWin([],[],fname);
imagesc(exWave,emWave,exemMatrix);
grid on;
identityLine;
xlabel('Excitation wave'); ylabel('Emission wave')
axis image

%%

Z = parafac(exemMatrix,1);
ieNewGraphWin; plot(emWave, Z{1}/max(Z{1}), 'k-','linewidth',1); grid on
xlabel('Wave (nm)'); ylabel('Emission energy')
title(sprintf('%s Parafac webfluor',fname));
% wavelength = emWave;
% data = Z{1};
% comment = 'Emission energy derived from parafac(webfluor datafile FAD.txt)';
% ieSaveSpectralFile(wavelength,data,comment,'webfluorFADEmission.mat');

ieNewGraphWin; plot(exWave,Z{2}/max(Z{2}), 'k-','linewidth',1); grid on
xlabel('Wave (nm)'); ylabel('Excitation sensitivity energy')
title(sprintf('%s Parafac webfluor',fname));
wavelength = exWave;
data = Z{2};
comment = 'Excitation energy derived from parafac(webfluor datafile FAD.txt)';
ieSaveSpectralFile(wavelength,data,comment,'webfluorFADExcitation.mat');


%% We may have multiple estimates of the same molecule
%
% These are scanned from the Monici 2005 paper
fName = fullfile(fiToolboxRootPath,'data','Monici','Elastin.mat');
f = fluorophoreRead(fName);
fluorophorePlot(f,'donaldson image');
%%
fluorophorePlot(f,'excitation');
fluorophorePlot(f,'emission');

%%

