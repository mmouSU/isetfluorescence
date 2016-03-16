function [ rootPath ] = fiToolboxRootPath()

rootPath=which('fiToolboxRootPath');
rootPath = fileparts(rootPath);

return
