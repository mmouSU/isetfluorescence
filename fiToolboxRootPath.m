function [ rootPath, parentPath ] = fiToolboxRootPath()

% [ rootPath ] = fiToolboxRootPath()
%
% Returns the absolute path to the root directory where the fiToolbox is
% installed.
%
% Copyright, Henryk Blasinski 2016

rootPath = which('fiToolboxRootPath');
rootPath = fileparts(rootPath);
parentPath = fileparts(rootPath);

return
