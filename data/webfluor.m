% Read data from the webfluor site
%
%   http://eemdb.uhnres.utoronto.ca/cgi-bin/WebObjects/WebFluor.woa/wo/jP5CODWCEJZC2BEKqPKzBn9aPzO/4.0.6.1.9.13.0.5.0.1.1#Flavin%20adenine%20dinucleotide%20[1]
%
% The site provides a bunch of Donaldson matrix data (which we are now
% calling Excitation-Emission Matrices).  You can write them out as a text
% file and then use this script to convert them to the isetfluorescence
% format.
%
% See also
%

%%
chdir(fullfile(fiToolboxRootPath,'local'));
fname = 'FAD.txt';

%% The  rows are excitation and the columns are emission
%
% The diagonal terms are the reflectance.  These are not quite zero.

T = readtable(fname);
exWave = T{:,1};
emWave = 260:5:750;
exemMatrix = T{:,2:end}';
mesh(exWave,emWave,exemMatrix)
xlabel('Excitation wave'); ylabel('Emission wave')

ieNewGraphWin;
imagesc(exWave,emWave,exemMatrix);
grid on;
identityLine;

%% Now save it out in the isetfluorescence data directory

