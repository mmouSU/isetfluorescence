% Read and analyze data from the webfluor site
%
% (This site went down in 2019).
%
% The data were measured by
%   Dacosta, R., Andersson, H., Wilson, B. "Molecular Fluorescence
%   Excitation?Emission Matrices Relevant to Tissue Spectroscopy"
%   Photochemistry and photobiology 2003/11/01, Vol 78, No 4. pp 384-392.
%
%   See pub for details about how the EEM data were "normalized" The graphs
%   in the pub plot Fluorescence Intensity (a.u.) - arbitrary units The
%   publication says that the data can be downloaded from
%   http://eemdb.uhnres.utoronto.ca/cgi-bin/WebObjects/WebFluor
%
%   **** By 10/05/2019, we could not connect to the website ***
%
%   http://eemdb.uhnres.utoronto.ca/cgi-bin/WebObjects/WebFluor.woa/wo/jP5CODWCEJZC2BEKqPKzBn9aPzO/4.0.6.1.9.13.0.5.0.1.1#Flavin%20adenine%20dinucleotide%20[1]
%
% http://eemdb.uhnres.utoronto.ca/cgi-bin/WebObjects/WebFluor.woa/wo/jP5uk9EI7aBc3ExSE6x1qodllrc/0.3
%
% The site provided Donaldson matrix data (which we are now calling
% Excitation-Emission Matrices).  You can write them out as a text file and
% then use this script to convert them to the isetfluorescence format.
%
% Not sure if the data are in units of energy or photons
%
% See also
%   parafac (isetfluor)
%

%%  Select a file from the webluor site

chdir(fullfile(fiToolboxRootPath,'data','sources','webfluor'));
txtFiles = dir('*.txt');

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

for ff = 1:numel(txtFiles)
    fname = txtFiles(ff).name;
    [~,fluorophoreName,e]= fileparts(fname);

    fprintf('Processing %s\n',fname);
    
    %% Read and plot the file
    %
    % Each column is the emission spectrum for a given excitation wavelength
    
    % The diagonal terms are the reflectance.  These are not quite zero.
    T = readtable(fname);
    ieNewGraphWin([],[],fname);
    exWave = T{:,1};       % Excitation wavelength range is in the 1st column
    emWave = 260:5:750;    % Emission wavelength range always this range
    
    % The excitation-emission matrix (Donaldson).  The first column is the
    % excitation wavelength samples.  Each row is the emission spectrum for
    % that wavelength.  The emission wavelength samples are 260:5:750,
    % according to the website (?) and paper (?).
    %
    % The transpose and flip of the table and exWave produces a matrix in which
    % the wavelengths increase across the columns (excitation) and down the
    % rows (emission).
    exemMatrix = fliplr(T{:,2:end}');
    exWave = flipud(exWave);
    mesh(exWave,emWave,exemMatrix);
    xlabel('Excitation wave'); ylabel('Emission wave')
    identityLine;
    title([fname,' original']);
    
    %% Make an image of the excitation-emission matrix
    
    %{
    ieNewGraphWin([],[],fname);
    imagesc(exWave,emWave,exemMatrix);
    identityLine;
    xlabel('Excitation wave'); ylabel('Emission wave')
    grid on; axis image
    %}
    
    %% Run parafac to estimate the excitation and emission spectra
    
    % We are using parafac to convert the Donaldson matrix into a single
    % fluorophore representation based only on excitation and emission spectra.
    % In general, parafac takes the excitation emission data from a substrate
    % and can interpret the data as being created from multiple fluorophores.
    % But in this case we think there is only a single fluorophore, and we use
    % parafac to obtain the ex and em vectors of a SINGLE fluorophore.
    %
    % We store these vectors in the isetfluor structure.
    %
    % A further simplification is that we require that the ex and em vectors be
    % represented on the same wavelength sampling points.  (Others would allow
    % them to be on different wavelength samples).
    %
    % When we are done, the excitation and emission spectra should be capable
    % of generating an excitation-emission matrix that approximates the
    % webfluor matrix.
    %
    
    % It is possible to specify the number of terms (fluorophores) represented
    % by a matrix.  We run this assuming 1.
    nFluorophores = 1;
    spectra = parafac(exemMatrix,nFluorophores);
    
    %%  Plot the excitation and emission spectra
    
    % JEF figured which output is which.  She compared the returns to published
    % data and found that they matched well. (December 2019).
    excitation = spectra{2};
    emission = spectra{1};
    
    %%Extrapolate to cover the same enlarged wavelength range
    
    % We extrapolate with zeros.
    minW = min(min(exWave),min(emWave));
    maxW = max(max(exWave),max(emWave));
    deltaW = 5;
    wave = minW:deltaW:maxW;
    emission   = interp1(emWave,emission,wave,'pchip',0);
    excitation = interp1(exWave,excitation,wave,'pchip',0);
    emission   = ieScale(emission,1);
    excitation = ieScale(excitation,1);
    
    ieNewGraphWin;
    plot(wave,emission,'r-',wave,excitation,'b-');
    grid on; title(fluorophoreName);
    
    %% Save out the data
    
    % We should check the solvent from the original paper
    thisF = fluorophoreCreate('type','custom',...
        'name',fluorophoreName,...
        'solvent','water', ...
        'wave', wave, ...
        'excitation',excitation,...
        'emission',emission);
    % disp(thisF)
    fluorophorePlot(thisF,'donaldson mesh')
    dMatrix = fluorophoreGet(thisF,'eem');
    
    %% Save out the fluorophore
    comment = sprintf('Webfluor data.  %s. Created %s',fname,datestr(now));
    saveName = fullfile(fiToolboxRootPath,'data','webfluor',fluorophoreName);
    fluorophoreSave(saveName,thisF,comment);
end

%% END


