function val = fluorescenceSceneGet(fl,param,varargin)
% Getter for the fluorescence structure
% 
% Copyright Henryk Blasinski, 2014

%% Parameter checking
if ~exist('fl','var') || isempty(fl), error('fluorescence structure required'); end
if ~exist('param','var') || isempty(param), error('param required'); end

val = [];

%% Main switch statement
param = ieParamFormat(param);
switch param
    case 'name'
        val = fl.name;
    case 'type'
        % Should always be 'illuminant'
        val = fl.type;
    case {'emission'}
        val = fl.data.emission;
        
    case {'excitation'}
        val = fl.data.excitation;
    case {'patchsize'}
        val = fl.data.patchSize;
        
    case {'illuminant'}
        val = fl.illuminant;
        
    case {'quantumefficiency'}
        val = fl.quantumEfficiency;
        
    case {'photons','Photons'}
                
        % Now we need to compute photons
        emissionCell = fluorescenceSceneGet(fl,'emission');
        excitationCell = fluorescenceSceneGet(fl,'excitation');
        [h, w] = size(emissionCell);
        
        wave = fluorescenceSceneGet(fl,'wave');
        deltaL = wave(2) - wave(1);
        nWaves = length(wave);
        
        emissionMat = zeros(h,w,nWaves);

        % Now we need to compute the excitation of each fluorophore
        ill = fluorescenceSceneGet(fl,'illuminant');
        ill = illuminantGet(ill,'photons');

        qe = fluorescenceSceneGet(fl,'quantumEfficiency');
        
        for x=1:w
            for y=1:h
                
                
                N = size(emissionCell{y,x},2);
                
                currEmission = emissionCell{y,x}*diag(1./sum(emissionCell{y,x}))/deltaL;
                currExcitation = excitationCell{y,x}*diag(1./max(excitationCell{y,x}));
                
                currEmission(isnan(currEmission)) = 0;
                currExcitation(isnan(currExcitation)) = 0;
                
                fullDonaldson = qe*currEmission*currExcitation'/N;
                donaldsonMat = tril(fullDonaldson,-1);
                
                
                % nAbsorbedPhotons = ill'*excitationCell{y,x}*deltaL;
                % nEmittedPhotons = sum(emissionCell{y,x})*deltaL;
                
                % scaleFactor = nAbsorbedPhotons./nEmittedPhotons;
                % scaleFactor(isnan(scaleFactor)) = 0;
                
                
                emissionMat(y,x,:) = donaldsonMat*ill*deltaL;
                
            end
        end

        patchSize = fluorescenceSceneGet(fl,'patchSize');
        photons = imageIncreaseImageRGBSize(emissionMat,patchSize);

        
        val = photons;
        
    case {'wave'}
        
        if isfield(fl,'spectrum'), val = fl.spectrum.wave;
        elseif ~isempty(varargin), val = sceneGet(varargin{1},'wave');
        end
        if isvector(val), val = val(:); end
        
    otherwise
        error('Unknown fluorescence scene parameter %s\n',param)
end

end
