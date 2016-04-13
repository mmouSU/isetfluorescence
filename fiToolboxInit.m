function fiToolboxInit()

close all;
clear all;
clc;

addpath(fullfile(fiToolboxRootPath));

addpath(fullfile(fiToolboxRootPath,'camera'));
addpath(fullfile(fiToolboxRootPath,'data'));

addpath(fullfile(fiToolboxRootPath,'estimation'));
addpath(fullfile(fiToolboxRootPath,'estimation','multistep'));
addpath(fullfile(fiToolboxRootPath,'estimation','nucnorm'));


addpath(fullfile(fiToolboxRootPath,'fluorescence'));
addpath(fullfile(fiToolboxRootPath,'io'));

addpath(fullfile(fiToolboxRootPath,'scripts'));
addpath(fullfile(fiToolboxRootPath,'scripts','validation'));
addpath(fullfile(fiToolboxRootPath,'scripts','simulations'));
addpath(fullfile(fiToolboxRootPath,'scripts','experiments'));
addpath(fullfile(fiToolboxRootPath,'scripts','plots'));

addpath(fullfile(fiToolboxRootPath,'utilities'));

end



