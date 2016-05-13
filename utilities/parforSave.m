function parforSave( fName, varargin )

% parforSave( fName, ... )
%
% A function used to save data from a parfor loop. This function DOES NOT
% work with Malab 2015b and newer. For example calling this function as
%    parforSave('test.mat',myVar);
% will save myVar in the test.mat file.
%
% Inputs:
%   fName - name of the .mat file where the variables will be saved 
%   
% Optional:
%   variables to be saved in the .mat file
%
% Copyright, The Internet


    for i=1:numel(varargin)
        eval([inputname(i+1) '=varargin{i};']);
    end
    
    save('-mat',fName,inputname(2));
    for i=2:numel(varargin)
        save('-mat',fName,inputname(i+1),'-append');
    end

end

