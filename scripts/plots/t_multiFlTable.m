% Generate a .tex file with a table contaning comparison results between
% multi-fluorophore and Suo et al. approaches. This script outputs
% Table 2 from the paper.
%
% Copyright, Henryk Blasinski 2016

close all;
clear all;
clc;

% Define the directory where figures will be saved. If saveDir = [], then
% figures are not saved.
% saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2');
saveDir = [];


nFilters = 8;
nChannels = 14;
nSamples = 24;

%% Load data
fName=  fullfile(fiToolboxRootPath,'results','experiments','multiFl_Macbeth+multiFl3.mat');
multiFl = load(fName);

fName=  fullfile(fiToolboxRootPath,'results','experiments','nucNorm_Macbeth+multiFl3.mat');
nucNorm = load(fName);


% Pixel values
measValsEst = multiFl.reflValsEst + multiFl.flValsEst;
[pixelError(1), pixelStd(1)] = fiComputeError(reshape(multiFl.measVals,[nFilters*nChannels nSamples]),reshape(measValsEst,[nFilters*nChannels nSamples]),'absolute');

measValsEst = nucNorm.reflValsEst + nucNorm.flValsEst;
[pixelError(2), pixelStd(2)] = fiComputeError(reshape(nucNorm.measVals,[nFilters*nChannels nSamples]),reshape(measValsEst,[nFilters*nChannels nSamples]),'absolute');

% Reflectance
[reflError(1), reflStd(1)] = fiComputeError(multiFl.reflEst,multiFl.reflRef,'absolute');

reflEst = [];
for i=1:nSamples
    reflEst = [reflEst diag(nucNorm.reflEst{i})];
end  
[reflError(2), reflStd(2)] = fiComputeError(reflEst,nucNorm.reflRef,'absolute');


% Donaldson matrix 
[dMatError(1), dMatStd(1)] = fiComputeError(multiFl.dMatEst,multiFl.dMatRef,'absolute');
[dMatError(2), dMatStd(2)] = fiComputeError(nucNorm.dMatEst,nucNorm.dMatRef,'absolute');

% Donaldson matrix 
[dMatScaledError(1), dMatScaledStd(1)] = fiComputeError(multiFl.dMatEst,multiFl.dMatRef,'normalized');
[dMatScaledError(2), dMatScaledStd(2)] = fiComputeError(nucNorm.dMatEst,nucNorm.dMatRef,'normalized');



%% Display results

fprintf('%20s | %20s | %20s | %20s | %20s \n','Algorithm','Pixels','Reflectance','Donaldson','Donaldson normalized');
fprintf('%s\n',repmat('-',[135 1]));
fprintf('%20s | %20f | %20f | %20f | %20s \n','Multi-fl.',pixelError(1),reflError(1),dMatError(1),dMatScaledError(1));
fprintf('%20s | %20f | %20f | %20f | %20f \n','Nuc. Norm',pixelError(2),reflError(2),dMatError(2),dMatScaledError(2));

%% Prepare a LaTeX table

str = '\\begin{table*}\n\\renewcommand{\\arraystretch}{1.3}\n\\centering\n';
str = [str '\\caption{Comparison between multi-fluorophore algorithms. The RMSE $\\pm$1 sd are shown. Suo et al. \\cite{Suo:14} represent reflectance as a square matrix.  Although their estimate sometimes includes non-zero off-diagonal terms, we estimate the RMSE reflectance error using only the diagonal terms.}\n'];
str = [str '\\label{tab:multiAccExp}\n\\begin{tabular}{| l | c | c | c | c |}\n\\hline\n' ...
            '\\multirow{2}{*}{\\diagbox{Algorithm}{Quantity}} & \\multirow{2}{*}{Pixel values} & '...
            '\\multirow{2}{*}{Reflectance} & \\multicolumn{2}{ c |}{Donaldson matrix} \\\\\n'];
str = [str '\\cline{4-5}\n& &  & {Absolute $\\times 10^{-2}$} & {Normalized} \\\\\n\\hline\n\\hline\n'];
str = [str sprintf('Ours - Multi-fluorophore & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ \\\\\\\\\n',pixelError(1),pixelStd(1),reflError(1),reflStd(1),dMatError(1)*100,dMatStd(1)*100,dMatScaledError(1),dMatScaledStd(1))];
str = [str sprintf('Suo et al. \\\\cite{Suo:14} & $%.2f \\\\pm %.2f$ & $ %.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $ %.2f\\\\pm %.2f$ \\\\\\\\\n',pixelError(2),pixelStd(2),reflError(2),reflStd(2),dMatError(2)*100,dMatStd(2)*100,dMatScaledError(2),dMatScaledStd(2))];   
str = [str '\\hline\n\\end{tabular}\n\\end{table*}\n'];

if ~isempty(saveDir)
    fName = fullfile(saveDir,'multiFlAccuracy.tex');
    fid = fopen(fName,'w');
    fprintf(fid,str);
    fclose(fid);
end



















