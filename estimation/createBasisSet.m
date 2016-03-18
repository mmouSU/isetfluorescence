function [ basis, score ] = createBasisSet( type, varargin )

p = inputParser;
p.addRequired('type',@ischar);
p.addParamValue('wave',400:10:700,@isnumeric);
p.addParamValue('n',5,@isscalar);

p.parse(type,varargin{:});
inputs = p.Results;

switch type
    case 'reflectance'
        
        fName = fullfile(isetRootPath,'data','surfaces','macbethChart');
        data = ieReadSpectra(fName,inputs.wave);

    case 'excitation'
        
        fName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell');
        fluorophores = fiReadFluorophoreSet(fName,'peakEmRange',[0 max(inputs.wave)],...
                                                  'peakExRange',[min(inputs.wave) Inf],...
                                                  'wave',inputs.wave);

        data = zeros(length(inputs.wave),length(fluorophores));
        for i=1:length(fluorophores)
            data(:,i) = fluorophoreGet(fluorophores(i),'normalized excitation');
        end

    case 'emission'
        
        fName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell');
        fluorophores = fiReadFluorophoreSet(fName,'peakEmRange',[0 max(inputs.wave)],...
                                                  'peakExRange',[min(inputs.wave) Inf],...
                                                  'wave',inputs.wave);

        data = zeros(length(inputs.wave),length(fluorophores));
        for i=1:length(fluorophores)
            data(:,i) = fluorophoreGet(fluorophores(i),'normalized emission');
        end   

end


[basis, ~, score] = pca(data','centered',false); 
basis = basis(:,1:inputs.n);
score = score(1:inputs.n);
end

