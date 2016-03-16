
ieInit;

data = load('./data/McNamaraBoswell.mat');

for i=1:length(data.McNamaraBoswell)
   
    sample = data.McNamaraBoswell(i);
    
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
    solName = regexprep(sample.solvent,'[^\d\w~!@#$%^&()_\-{}.]*',''); 
    
    fName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell',[cName '+' solName '.mat']);
    
        
    fiSaveFluorophore(fName,fl,'McNamara-Boswell dataset');
    
    
end


