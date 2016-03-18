function [ val ] = fluorescentSceneGet(flScene,param,varargin)

%% Parameter checking
if ~exist('flScene','var') || isempty(flScene), error('Fluorescent scene structure required'); end
if ~exist('param','var') || isempty(param), error('param required'); end

val = [];

%% Main switch statement
param = lower(param);
param = strrep(param,' ','');

switch param
    case 'name'
        val = flScene.name;

    case 'type'
        % Should always be 'fluorescent scene'
        val = flScene.type;

    case 'nfluorophores'
        if ~checkfields(flScene,'fluorophores'), val = []; return; end
        val = size(flScene.fluorophores,3);

    case 'nwaves'
        val = length(flScene.spectrum.wave);
       
    case {'donaldsonmatrix'}
        
        nFluorophores = fluorescentSceneGet(flScene,'nFluorophores');

        if length(varargin) == 3 % we've passed a vector of wavelengths
            wave = varargin{3};
        else
            wave = fluorescentSceneGet(flScene,'wave');
        end

        nWave = length(wave);
        DM = zeros(nWave);
        for i=1:nFluorophores

            fl = fluorophoreSet(flScene.fluorophores(varargin{1},varargin{2},i),'wave',wave);
            DM = DM + fluorophoreGet(fl,'Donaldson matrix');
        end

        val = DM;

    case 'scenesize'
        val = [flScene.height; flScene.width];

    case 'size'
        val = [flScene.rows; flScene.cols];


    case {'photons'}
        illWave  = illuminantGet(varargin{1},'wave');
        nWaves = length(illWave);
        illSpd = illuminantGet(varargin{1},'photons');
        

        photons = zeros(flScene.rows,flScene.cols,nWaves);

        for xx=1:flScene.cols
            for yy=1:flScene.rows
                DM = fluorescentSceneGet(flScene,'Donaldson matrix',yy,xx,illWave);
                photons(yy,xx,:) = DM*illSpd;
            end
        end

        val = imresize(photons,[flScene.height,flScene.width],'nearest');
        
    case {'wave'}
        val = flScene.spectrum.wave;

    case {'donaldsonreference'}
        
        sz = fluorescentSceneGet(flScene,'size');
        val = cell(sz(1)*sz(2),1);
        for xx=1:sz(2)
            for yy=1:sz(1)
                DM = fluorescentSceneGet(flScene,'Donaldson matrix',yy,xx);
                val{(xx-1)*sz(1) + yy} = DM;
            end
        end
    
    case {'excitationreference'}
        sz = fluorescentSceneGet(flScene,'size');
        nFluorophores = fluorescentSceneGet(flScene,'nFluorophores');
        nWaves = fluorescentSceneGet(flScene,'nWaves');


        val = zeros(nWaves,sz(1)*sz(2)*nFluorophores);
        for xx=1:sz(2)
            for yy=1:sz(1)
                for zz=1:nFluorophores
                    val(:,(xx-1)*nFluorophores*sz(1) + (yy-1)*nFluorophores + zz) = fluorophoreGet(flScene.fluorophores(yy,xx,zz),'excitation');
                end
            end
        end 

    case {'emissionreference'}
        sz = fluorescentSceneGet(flScene,'size');
        nFluorophores = fluorescentSceneGet(flScene,'nFluorophores');
        nWaves = fluorescentSceneGet(flScene,'nWaves');


        val = zeros(nWaves,sz(1)*sz(2)*nFluorophores);
        for xx=1:sz(2)
            for yy=1:sz(1)
                for zz=1:nFluorophores
                    qe = fluorophoreGet(flScene.fluorophores(yy,xx,zz),'qe');
                    deltaL = fluorophoreGet(flScene.fluorophores(yy,xx,zz),'deltaWave');
                    val(:,(xx-1)*nFluorophores*sz(1) + (yy-1)*nFluorophores + zz) = deltaL*qe*fluorophoreGet(flScene.fluorophores(yy,xx,zz),'emission');
                end
            end
        end 

    otherwise
        error('Unknown fluorescent scene parameter %s\n',param)
end

end

