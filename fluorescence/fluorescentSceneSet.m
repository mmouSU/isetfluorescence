function [ flScene ] = fluorescentSceneSet( flScene, param, val, varargin)
%%

if ~exist('flScene','var') || isempty(flScene), error('Fluorophore structure required'); end
if ~exist('param','var') || isempty(param), error('param required'); end
if ~exist('val','var') , error('val is required'); end

%%
param = lower(param);
param = strrep(param,' ','');

switch param
    case 'name'
        flScene.name = val;
        
    case 'type'
        if ~strcmpi(val,'fluorescent scene'), error('Type must be ''fluorescent scene'''); end
        flScene.type = val;

    case 'scenesize'
        flScene.height = val(1);
        flScene.width = val(2);

    case 'size'
        flScene.rows = val(1);
        flScene.cols = val(2);
        
    case 'qe'
        % We 'split' the quantum efficiency evenly across all teh fluorophores at a particular spatial location. 
        if (val > 1), warning('Qe greater than one, truncating to 1'); end
        val = min(max(val,0),1);
        
        allFluorophores = length(flScene.fluorophores(:));
        nFluorophores = size(flScene.fluorophores,3);

        
        for i=1:allFluorophores
            if isscalar(val)
                % If we provide a scalr qe then we set the same value at
                % all spatial locations
                flScene.fluorophores(i) = fluorophoreSet(flScene.fluorophores(i),'qe',val/nFluorophores);
            else
                % If it is a matrix, then each fluorophore get set a
                % specific value
                flScene.fluorophores(i) = fluorophoreSet(flScene.fluorophores(i),'qe',val(i));
            end
        end
        
    case 'fluorophores'
        flScene.fluorophores = val;
      
        
    otherwise
        error('Unknown fluorophore parameter %s\n',param)
end


end

