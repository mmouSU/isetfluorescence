% Generate a .tex file with a table contaning comparison results between
% multi-fluorophore and Suo et al. approaches on simulated data from the
% McNamara-Boswell fluorophore dataset. This script outputs
% Table 1 from the supplemental material.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
saveDir = fullfile('~','Dropbox','Fluorescence Imaging','Notes','TPAMI - 3rd attempt');
% saveDir = [];


nFilters = 8;
nChannels = 14;
nSamples = 24;

%% Load data
fName = fullfile(fiToolboxRootPath,'results','evaluation','McNamara-Boswell_simCompare_multiFl.mat');
load(fName);


% Pixel values
pixelError(1) = mean(multiFlTotalPixelErr);
pixelStd(1) = mean(multiFlTotalPixelStd);

% The nuclear norm algorithm did not always converge. We assume a somewhat
% arbitrary detection threshold.
valid = nucNormTotalPixelErr <= 0.5;
pixelError(2) = mean(nucNormTotalPixelErr(valid));
pixelStd(2) = mean(nucNormTotalPixelStd(valid));


% Reflectance
reflError(1) = mean(multiFlReflErr);
reflStd(1) = mean(multiFlReflStd);

reflError(2) = mean(nucNormReflErr(valid));
reflStd(2) = mean(nucNormReflStd(valid));


% Donaldson matrix (absolute) 

dMatError(1) = mean(multiFldMatErr);
dMatStd(1) = mean(multiFldMatStd);

dMatError(2) = mean(nucNormdMatErr(valid));
dMatStd(2) = mean(nucNormdMatStd(valid));

% Donaldson matrix (normalized)
dMatScaledError(1) = mean(multiFldMatNormErr);
dMatScaledStd(1) = mean(multiFldMatNormStd);

dMatScaledError(2) = mean(nucNormdMatNormErr(valid));
dMatScaledStd(2) = mean(nucNormdMatNormStd(valid));



%% Display results

fprintf('%20s | %20s | %20s | %20s | %20s \n','Algorithm','Pixels','Reflectance','Donaldson','Donaldson normalized');
fprintf('%s\n',repmat('-',[135 1]));
fprintf('%20s | %20f | %20f | %20f | %20s \n','Multi-fl.',pixelError(1),reflError(1),dMatError(1),dMatScaledError(1));
fprintf('%20s | %20f | %20f | %20f | %20f \n','Nuc. Norm',pixelError(2),reflError(2),dMatError(2),dMatScaledError(2));

%% Prepare a LaTeX table

str = '\\begin{table*}\n\\renewcommand{\\arraystretch}{1.3}\n\\centering\n';
str = [str '\\caption{Comparison between multi-fluorophore algorithms on simulated data from the McNamara-Boswell dataset. The RMSE $\\pm$1 sd are shown. Suo et al. \\cite{Suo:14} represent reflectance as a square matrix.  Although their estimate sometimes includes non-zero off-diagonal terms, we estimate the RMSE reflectance error using only the diagonal terms.}\n'];
str = [str '\\label{tab:multiAccExp}\n\\begin{tabular}{| l | c | c | c | c |}\n\\hline\n' ...
            '\\multirow{2}{*}{\\diagbox{Algorithm}{Quantity}} & \\multirow{2}{*}{Pixel values} & '...
            '\\multirow{2}{*}{Reflectance} & \\multicolumn{2}{ c |}{Donaldson matrix} \\\\\n'];
str = [str '\\cline{4-5}\n& &  & {Absolute $\\times 10^{-2}$} & {Normalized} \\\\\n\\hline\n\\hline\n'];
str = [str sprintf('Ours - Multi-fluorophore & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ \\\\\\\\\n',pixelError(1),pixelStd(1),reflError(1),reflStd(1),dMatError(1)*100,dMatStd(1)*100,dMatScaledError(1),dMatScaledStd(1))];
str = [str sprintf('Suo et al. \\\\cite{Suo:14} & $%.2f \\\\pm %.2f$ & $ %.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $ %.2f\\\\pm %.2f$ \\\\\\\\\n',pixelError(2),pixelStd(2),reflError(2),reflStd(2),dMatError(2)*100,dMatStd(2)*100,dMatScaledError(2),dMatScaledStd(2))];   
str = [str '\\hline\n\\end{tabular}\n\\end{table*}\n'];

if ~isempty(saveDir)
    fName = fullfile(saveDir,'multiFlSimAccuracy.tex');
    fid = fopen(fName,'w');
    fprintf(fid,str);
    fclose(fid);
end



















