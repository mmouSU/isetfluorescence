% Generate a .tex file with a table contaning comparison results between
% CIM, single fluorophore and the multistep appraches. This script outputs
% Table 3 from the paper.
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
fName=  fullfile(fiToolboxRootPath,'results','experiments','em_Macbeth+fl.mat');
em = load(fName);

fName=  fullfile(fiToolboxRootPath,'results','experiments','fl_Macbeth+fl.mat');
fl = load(fName);

fName=  fullfile(fiToolboxRootPath,'results','experiments','multistep_Macbeth+fl.mat');
multi = load(fName);


% Pixel values
measValsEst = em.reflValsEst + em.flValsEst;
[pixelError(1), pixelStd(1)] = fiComputeError(reshape(em.measVals,[nFilters*nChannels nSamples]),reshape(measValsEst,[nFilters*nChannels nSamples]),'absolute');

measValsEst = fl.reflValsEst + fl.flValsEst;
[pixelError(2), pixelStd(2)] = fiComputeError(reshape(fl.measVals,[nFilters*nChannels nSamples]),reshape(measValsEst,[nFilters*nChannels nSamples]),'absolute');

measValsEst = multi.reflValsEst + multi.flValsEst;
[pixelError(3), pixelStd(3)] = fiComputeError(reshape(multi.measVals,[nFilters*nChannels nSamples]),reshape(measValsEst,[nFilters*nChannels nSamples]),'absolute');


% Reflectance
[reflError(1), reflStd(1)] = fiComputeError(em.reflEst,em.reflRef,'absolute');
[reflError(2), reflStd(2)] = fiComputeError(fl.reflEst,fl.reflRef,'absolute');
[reflError(3), reflStd(3)] = fiComputeError(multi.reflEst,multi.reflRef,'absolute');


% Emission (normalized)
[emError(1), emStd(1)] = fiComputeError(em.emEst,em.emRef,'normalized');
[emError(2), emStd(2)] = fiComputeError(fl.emEst,fl.emRef,'normalized');
[emError(3), emStd(3)] = fiComputeError(multi.emEst,multi.emRef,'normalized');


% Excitation (normalized)
[exError(2), exStd(2)] = fiComputeError(fl.exEst,fl.exRef,'normalized');
[exError(3), exStd(3)] = fiComputeError(multi.exEst,multi.exRef,'normalized');

% Excitation (scaled)
[exScaledError(2), exScaledStd(2)] = fiComputeError(fl.exEst*diag(max(fl.emEst)),fl.exRef*diag(max(fl.emRef)),'absolute');
[exScaledError(3), exScaledStd(3)] = fiComputeError(multi.exEst*diag(max(multi.emEst)),multi.exRef*diag(max(multi.emRef)),'absolute');


%% Display results

fprintf('%20s | %20s | %20s | %20s | %20s | %20s \n','Algorithm','Pixels','Reflectance','Emission','Excitation','Excitaiton scaled');
fprintf('%s\n',repmat('-',[135 1]));
fprintf('%20s | %20f | %20f | %20f | %20s | %20s \n','Emission',pixelError(1),reflError(1),emError(1),'-','-');
fprintf('%20s | %20f | %20f | %20f | %20f | %20f \n','Single fl.',pixelError(2),reflError(2),emError(2),exError(2),exScaledError(2));
fprintf('%20s | %20f | %20f | %20f | %20f | %20f \n','Multistep',pixelError(3),reflError(3),emError(3),exError(3),exScaledError(3));

%% Prepare a LaTeX table

str = '\\begin{table*}\n\\renewcommand{\\arraystretch}{1.3}\n\\centering\n';
str = [str '\\caption{Comparison of single fluorophore estimation algorithms. The RMSE $\\pm$ 1 sd are shown.'...
    'There is a scaling ambiguity between the excitation and emission spectra; without loss of generality we assume' ...
    ' that the peak emission spectrum is scaled to a value of 1 and that the fluorescence intensity is contained entirely '...
    'within the scale of the excitation spectrum. Hence, we include the RMSE of the excitation estimate with and without scaling, '...
    'but there is no need to do so for the emission spectrum.}\n'];
str = [str '\\label{tab:singleAccExp}\n\\begin{tabular}{| l | c | c | c | c | c |}\n\\hline\n' ...
            '\\multirow{2}{*}{\\diagbox{Algorithm}{Quantity}} & \\multirow{2}{*}{Pixel values} & '...
            '\\multirow{2}{*}{Reflectance} & \\multirow{2}{*}{Emission} & \\multicolumn{2}{ c |}{Excitation} \\\\\n'];
str = [str '\\cline{5-6}\n& &  &  & {Scale $\\times 10^{-1}$} & {Shape} \\\\\n\\hline\n\\hline\n'];
str = [str sprintf('Ours - CIM & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & -- & -- \\\\\\\\\n',pixelError(1),pixelStd(1),reflError(1),reflStd(1),emError(1),emStd(1))];
str = [str sprintf('Ours - Single fluorophore & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ \\\\\\\\\n',pixelError(2),pixelStd(2),reflError(2),reflStd(2),emError(2),emStd(2),exScaledError(2)*10,exScaledStd(2)*10,exError(2),exStd(2))];
str = [str sprintf('Fu et al. \\\\cite{Fu:14} & $%.2f \\\\pm %.2f$ & $ %.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $%.2f \\\\pm %.2f$ & $ %.2f\\\\pm %.2f$ \\\\\\\\\n',pixelError(3),pixelStd(3),reflError(3),reflStd(3),emError(3),emStd(3),exScaledError(3)*10,exScaledStd(3)*10,exError(3),exStd(3))];    
str = [str '\\hline\n\\end{tabular}\n\\end{table*}\n'];

if ~isempty(saveDir)
    fName = fullfile(saveDir,'singleFlAccuracy.tex');
    fid = fopen(fName,'w');
    fprintf(fid,str);
    fclose(fid);
end



















