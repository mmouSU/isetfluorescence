% Read data from the webfluor site
%
%   http://eemdb.uhnres.utoronto.ca/cgi-bin/WebObjects/WebFluor.woa/wo/jP5CODWCEJZC2BEKqPKzBn9aPzO/4.0.6.1.9.13.0.5.0.1.1#Flavin%20adenine%20dinucleotide%20[1]
%
% http://eemdb.uhnres.utoronto.ca/cgi-bin/WebObjects/WebFluor.woa/wo/jP5uk9EI7aBc3ExSE6x1qodllrc/0.3
%
% The site provides a bunch of Donaldson matrix data (which we are now
% calling Excitation-Emission Matrices).  You can write them out as a text
% file and then use this script to convert them to the isetfluorescence
% format.
%
% Not sure if the data are in units of energy or photons
%
% See also
%

%%
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));
%{
fname = 'FAD.txt';
fname = 'Flavin.txt;
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


%% The  rows are excitation and the columns are emission
%
% The diagonal terms are the reflectance.  These are not quite zero.

T = readtable(fname);
exWave = T{:,1};
emWave = 260:5:750;
exemMatrix = T{:,2:end}';
mesh(exWave,emWave,exemMatrix)
xlabel('Excitation wave'); ylabel('Emission wave')
identityLine
title(fname);

%%
ieNewGraphWin([],[],fname);
imagesc(exWave,emWave,exemMatrix);
grid on;
identityLine;
xlabel('Excitation wave'); ylabel('Emission wave')
axis image

%% Now save it out in the isetfluorescence data directory

