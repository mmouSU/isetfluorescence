function [ val ] = fluorescentSceneGet(flScene,param,varargin)

p = inputParser;
p.addRequired('flScene');
p.addRequired('param',@ischar);
p.addParamValue('sceneSize',[0 0],@isvector);
p.addParamValue('flPos',[1 1],@isvector);
p.addParamValue('wave',[],@isnumeric);
p.addParamValue('illuminant',[],@isstruct);

p.parse(flScene,param,varargin{:});
inputs = p.Results;


val = [];

%% Main switch statement
inputs.param = lower(inputs.param);
inputs.param = strrep(inputs.param,' ','');

switch inputs.param
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

        if isempty(inputs.wave)  % we've passed a vector of wavelengths
            inputs.wave = fluorescentSceneGet(flScene,'wave');
        end

        nWave = length(inputs.wave);
        DM = zeros(nWave);
        for i=1:nFluorophores

            fl = fluorophoreSet(flScene.fluorophores(inputs.flPos(1),inputs.flPos(2),i),'wave',inputs.wave);
            DM = DM + fluorophoreGet(fl,'Donaldson matrix');
        end

        val = DM;


    case 'size'
        val = [size(flScene.fluorophores,1); size(flScene.fluorophores,2)];


    case {'photons'}
        illWave  = illuminantGet(inputs.illuminant,'wave');
        nWaves = length(illWave);
        illSpd = illuminantGet(inputs.illuminant,'photons');
        
        sz = fluorescentSceneGet(flScene,'size');

        val = zeros(sz(1),sz(2),nWaves);

        for xx=1:sz(2)
            for yy=1:sz(1)
                DM = fluorescentSceneGet(flScene,'Donaldson matrix','flPos',[yy,xx],'wave',illWave);
                val(yy,xx,:) = DM*illSpd;
            end
        end

        
    case {'wave'}
        val = flScene.spectrum.wave;

    case {'donaldsonreference'}
        
        % For all reference calls we can choose the scene size, so that if
        % say we have a macbeth chart an a uniform fluorophore we will get
        % reference data for each Macbeht patch.

        sz = fluorescentSceneGet(flScene,'size');
        if sum(inputs.sceneSize) == 0
            sceneSize = sz;
        else
            sceneSize = inputs.sceneSize;
        end
        [mapX, mapY] = mapToMatchSceneSize(sz(1),sz(2),sceneSize(1),sceneSize(2));


        val = cell(sceneSize(1)*sceneSize(2),1);
        for xx=1:sceneSize(2)
            for yy=1:sceneSize(1)
                DM = fluorescentSceneGet(flScene,'Donaldson matrix','flPos',[mapY(yy,xx), mapX(yy,xx)]);
                val{(xx-1)*sceneSize(1) + yy} = DM;
            end
        end
    
    case {'excitationreference'}
        sz = fluorescentSceneGet(flScene,'size');
        if sum(inputs.sceneSize) == 0
            sceneSize = sz;
        else
            sceneSize = inputs.sceneSize;
        end
        [mapX, mapY] = mapToMatchSceneSize(sz(1),sz(2),sceneSize(1),sceneSize(2));



        nFluorophores = fluorescentSceneGet(flScene,'nFluorophores');
        nWaves = fluorescentSceneGet(flScene,'nWaves');


        val = zeros(nWaves,sceneSize(1)*sceneSize(2)*nFluorophores);
        for xx=1:sceneSize(2)
            for yy=1:sceneSize(1)

                locX = mapX(yy,xx);
                locY = mapY(yy,xx);

                for zz=1:nFluorophores
                    val(:,(xx-1)*nFluorophores*sceneSize(1) + (yy-1)*nFluorophores + zz) = fluorophoreGet(flScene.fluorophores(locY,locX,zz),'excitation');
                end
            end
        end 

    case {'emissionreference'}
        sz = fluorescentSceneGet(flScene,'size');
        if sum(inputs.sceneSize) == 0
            sceneSize = sz;
        else
            sceneSize = inputs.sceneSize;
        end
        [mapX, mapY] = mapToMatchSceneSize(sz(1),sz(2),sceneSize(1),sceneSize(2));
        nFluorophores = fluorescentSceneGet(flScene,'nFluorophores');
        nWaves = fluorescentSceneGet(flScene,'nWaves');


        val = zeros(nWaves,sceneSize(1)*sceneSize(2)*nFluorophores);
        for xx=1:sceneSize(2)
            for yy=1:sceneSize(1)

                locX = mapX(yy,xx);
                locY = mapY(yy,xx);

                for zz=1:nFluorophores
                    qe = fluorophoreGet(flScene.fluorophores(locY,locX,zz),'qe');
                    deltaL = fluorophoreGet(flScene.fluorophores(locY,locX,zz),'deltaWave');
                    val(:,(xx-1)*nFluorophores*sceneSize(1) + (yy-1)*nFluorophores + zz) = deltaL*qe*fluorophoreGet(flScene.fluorophores(locY,locX,zz),'emission');
                end
            end
        end 

    otherwise
        error('Unknown fluorescent scene parameter %s\n',param)
end

end


function [mapX, mapY] = mapToMatchSceneSize(flHeight,flWidth,scHeight,scWidth)

if flHeight > scHeight || flWidth > scWidth, error('Fluorescent scene has more test patches than the input scene'); end

if mod(scHeight,flHeight) ~= 0 || mod(scWidth,flWidth) ~=0, error('Fluorescent and reflective scene patches are not multiples of one another'); end

deltaX = scWidth/flWidth;
deltaY = scHeight/flHeight;

mapX = zeros(scHeight,scWidth);
mapY = zeros(scHeight,scWidth);

for x=1:flWidth
    for y=1:flHeight

        mapX((y-1)*deltaY+1:y*deltaY,(x-1)*deltaX+1:x*deltaX) = x;
        mapY((y-1)*deltaY+1:y*deltaY,(x-1)*deltaX+1:x*deltaX) = y;

    end
end


end

