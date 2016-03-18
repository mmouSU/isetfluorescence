
ieInit;

data = load('./data/LifeTechnologies.mat');

for i=1:length(data.LifeTechnologies)
   
    sample = data.LifeTechnologies(i);
    fl = [];
    if ~isempty(sample.excitation)
        fl = fluorophoreCreate('type','custom',...             
                                'wave',sample.wave,...
                                'name',sample.dye,...
                                'solvent',sample.solvent,...
                                'excitation',sample.excitation,...
                                'emission',sample.emission);
    else
        fl = fluorophoreCreate('type','custom',...             
                            'wave',sample.wave,...
                            'name',sample.dye,...
                            'solvent',sample.solvent,...
                            'excitation',sample.absorption,...
                            'emission',sample.emission);
    end
                        
    cName = regexprep(sample.dye,'[^\d\w~!@#$%^&()_\-{}.]*','');                  
    % solName = regexprep(sample.solvent,'[^\d\w~!@#$%^&()_\-{}.]*',''); 
    
    % fName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell',[cName '+' solName '.mat']);
    fName = fullfile(fiToolboxRootPath,'data','LifeTechnologies',[cName '.mat']);
        
    fiSaveFluorophore(fName,fl,'Life Technologies dataset');
    
    
end


