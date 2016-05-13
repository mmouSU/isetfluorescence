function [ rootPath ] = fiToolboxRootPath()

% [ rootPath ] = fiToolboxRootPath()
%
% Returns the absolute path to the root directory where the fiToolbox is
% installed.
%
% Copyright, Henryk Blasinski 2016

rootPath = which('fiToolboxRootPath');
rootPath = fileparts(rootPath);

return
