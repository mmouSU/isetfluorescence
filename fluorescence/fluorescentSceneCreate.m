function [ flScene, fluorophoreIDs ] = fluorescentSceneCreate( varargin )

p = inputParser;

p.addParamValue('type','Default',@ischar);
p.addParamValue('wave',400:10:700,@isvector);
p.addParamValue('name','Default',@ischar);
p.addParamValue('height',1,@isscalar);
p.addParamValue('width',1,@isscalar);
p.addParamValue('stokesShiftRange',[0 Inf],@isvector);
p.addParamValue('peakEmRange',[420 680],@isvector);
p.addParamValue('peakExRange',[420 680],@isvector);
p.addParamValue('nFluorophores',1,@isscalar);
p.addParamValue('qe',1,@isscalar);

p.parse(varargin{:});
inputs = p.Results;

flScene.type = 'fluorescent scene';
flScene.name = inputs.name;
flScene = initDefaultSpectrum(flScene,'custom',inputs.wave);

switch inputs.type
    
    otherwise
        % Create a default fluorescent scene
        setName = fullfile(fiToolboxRootPath,'data','LifeTechnologies');
        [flSet, fluorophoreIDs] = fiReadFluorophoreSet(setName,'wave',inputs.wave,...
            'stokesShiftRange',inputs.stokesShiftRange,...
            'peakEmRange',inputs.peakEmRange,...
            'peakExRange',inputs.peakExRange);
        
        nFluorophores = inputs.height*inputs.width*inputs.nFluorophores;
        
        ids = randi(length(flSet),nFluorophores,1);
        selFl = flSet(ids);
        fluorophoreIDs = fluorophoreIDs(ids);
        
        flScene = fluorescentSceneSet(flScene,'fluorophores',reshape(selFl,[inputs.height, inputs.width, inputs.nFluorophores]));        
        flScene = fluorescentSceneSet(flScene,'qe',inputs.qe/inputs.nFluorophores);
        flScene = fluorescentSceneSet(flScene,'size', [inputs.height inputs.width]);
        flScene = fluorescentSceneSet(flScene,'scene size',[inputs.height inputs.width]);
        
end



end

