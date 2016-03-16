function [ fl ] = fiReadFluorophore( fName, varargin )

p = inputParser;
p.addRequired('fName',@ischar);
p.addParamValue('wave','',@isvector);
p.addParamValue('qe',1,@isscalar);

p.parse(fName,varargin{:});
inputs = p.Results;


data = load(inputs.fName);

fl = fluorophoreCreate('type','custom',...
                       'wave',data.wave,...
                       'name',data.name,...
                       'solvent',data.solvent,...
                       'excitation',data.excitation,...
                       'emission',data.emission);
                   
fl = fluorophoreSet(fl,'qe',inputs.qe);
if ~isempty(inputs.wave)
    fl = fluorophoreSet(fl,'wave',inputs.wave);
end



end

