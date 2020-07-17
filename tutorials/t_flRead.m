%% Here is the way to create a Donaldson matrix from fluorescent data
%
% The data are in the isetfluorescence repository.  They are stored as a
% fluorophore structure.  There are many in the data/ directories.
%

%% Write out the text file

% Set the wavelength sampling you would like
wave = 395:10:705; 

% Here is one fluorophore, read in with that wavelength sampling.
FAD = fluorophoreRead('FAD','wave',wave);

% Here is the excitation emission matrix
eem = fluorophoreGet(FAD,'eem');
%{
 fluorophorePlot(FAD,'donaldson mesh');
%}
%{
 dWave = fluorophoreGet(FAD,'delta wave');
 wave = fluorophoreGet(FAD,'wave');
 ex = fluorophoreGet(FAD,'excitation');
 ieNewGraphWin; 
 plot(wave,sum(eem)/dWave,'k--',wave,ex/max(ex(:)),'r:')
%}

% The data are converted to a vector like this
wave = fluorophoreGet(FAD,'wave');
flatEEM = eem';
vec = [wave(1) wave(2)-wave(1) wave(end) flatEEM(:)'];

%% To write it out the values as a txt file vector

fid = fopen('testFluorophore.txt','w');
fprintf(fid,'%f ',vec);
fclose(fid);

%% Check the code

fid = fopen('testFluorophore.txt','r');
txt = textscan(fid,'%f');
fclose(fid);

val = txt{1};
wList = val(1:3);
wave = wList(1):wList(2):wList(3);
eem = val(4:end);
eem = reshape(eem,32,32);

ieNewGraphWin; imagesc(wave,wave,eem');
axis image

%% END