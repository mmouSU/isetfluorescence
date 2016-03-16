function [ flSet, indx ] = fiReadFluorophoreSet( dirName, varargin )

p = inputParser;
p.addRequired('dirName',@ischar);
p.addParamValue('stokesShiftRange',[0 Inf],@isvector);
p.addParamValue('peakEmRange',[0 Inf],@isvector);
p.addParamValue('peakExRange',[0 Inf],@isvector);
p.addParamValue('wave','',@isvector);
p.addParamValue('verbose',false,@islogical);

p.parse(dirName,varargin{:});
inputs = p.Results;


files = dir(fullfile(inputs.dirName,'*.mat'));
nFiles = length(files);

flSet = [];
indx = [];
for i=1:nFiles
    
    fName = fullfile(dirName,files(i).name);
    fl = fiReadFluorophore(fName);
    
    % Remove the fluorphores outside the Stokes shift range
    if fluorophoreGet(fl,'Stokes shift') < min(inputs.stokesShiftRange) ||  fluorophoreGet(fl,'Stokes shift') > max(inputs.stokesShiftRange)
        continue;
    end
    
    % Remove the fluorphores outside the emission range
    if fluorophoreGet(fl,'Peak emission') < min(inputs.peakEmRange) ||  fluorophoreGet(fl,'Peak emission') > max(inputs.peakEmRange)
        continue;
    end
    
    % Remove the fluorphores outside the emission range
    if fluorophoreGet(fl,'Peak excitation') < min(inputs.peakExRange) ||  fluorophoreGet(fl,'Peak excitation') > max(inputs.peakExRange)
        continue;
    end
    
    if ~isempty(inputs.wave)
        fl = fluorophoreSet(fl,'wave',inputs.wave);
    end
    
    if inputs.verbose
        fprintf('%s : peak excitation %i, peak emission %i, Stokes shift %i\n',fluorophoreGet(fl,'name'),...
                                                                               fluorophoreGet(fl,'Peak excitation'),...
                                                                               fluorophoreGet(fl,'Peak emission'),...
                                                                               fluorophoreGet(fl,'Stokes shift'));
    end
    
    flSet = [flSet; fl];
    indx = [indx; i];
    
end


end

