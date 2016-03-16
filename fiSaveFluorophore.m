function fiSaveFluorophore( fName, fl, comment )

p = inputParser;
p.addRequired('fName',@ischar);
p.addRequired('fl',@isstruct);
p.addOptional('comment','',@ischar);

p.parse(fName,fl,comment);
inputs = p.Results;

name = inputs.fl.name;
solvent = inputs.fl.solvent;
excitation = inputs.fl.excitation/max(inputs.fl.excitation);
emission = inputs.fl.emission/max(inputs.fl.emission);
comment = inputs.comment;
wave = inputs.fl.spectrum.wave;

save(fName,'name','solvent','excitation','emission','comment','wave');

end

