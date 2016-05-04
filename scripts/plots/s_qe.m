close all;
clear all;
clc;

fName = fullfile(fiToolboxRootPath,'results','evaluation','McNamara-Boswell_simQe_Fl.mat');
load(fName);

figure;
hold on; grid on; box on;
errorbar(flQe, totalPixelErr,totalPixelStd/sqrt(24),'-d');

errorbar(flQe, reflErr,reflStd/sqrt(24));
errorbar(flQe, emNormErr,emNormStd/sqrt(24));
errorbar(flQe, exNormErr,exNormStd/sqrt(24));

set(gca,'xscale','log');
set(gca,'yscale','log');

