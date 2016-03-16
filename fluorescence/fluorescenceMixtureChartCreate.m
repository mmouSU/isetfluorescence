function [fChart, absorption, emission, weights, names, IDs] = fluorescenceMixtureChartCreate(patchSize,h,w,spectrum,ill,stokesLimit,nCompounds, compoundIDs, weights, qe)

fChart = fluorescenceSceneCreate('Fluorescence chart',spectrum);

% fluorescenceChartObject.name = 'Fluorescence chart';
% fluorescenceChartObject.type = 'fluorescence scene';

% This is the size in pixels of each testchart patch
if ieNotDefined('patchSize'), patchSize = 16;   end

% These are the limits on the Stokes shifts of the compounds which will
% form the test chart
if ieNotDefined('stokesLimit'), stokesLimit = [0 100]; end

% This parameter specifies the quantuum efficiency, i.e. the ratio between
% the number of emited to absorbed photons.
if ieNotDefined('qe'), qe = 1; end

% Define the size of the chart
if ieNotDefined('h'), h = 4; end
if ieNotDefined('w'), w = 6; end

if ieNotDefined('nCompounds'), nCompounds = 10; end
if ieNotDefined('qe'), qe = ones(nCompounds,1); end


deltaL = spectrum(2) - spectrum(1);
nWaves = length(spectrum);
nPatches = h*w;
fName = fullfile(isetRootPath,'data','surfaces','McNamaraBoswell.mat');
load(fName);

% If we don't specify patches, we draw them at random from the databse
if ieNotDefined('compoundIDs')
        
    IDs = 1:length(McNamaraBoswell);
    stokes = ([McNamaraBoswell(:).stokesShift] >= stokesLimit(1)) & ([McNamaraBoswell(:).stokesShift] < stokesLimit(2));
    
    % Restrict the absortpion peak
    t1 = {McNamaraBoswell(:).maxAbs};
    t1indx = cellfun(@isempty,t1);
    t1(t1indx) = {0};
    maxAbs = cell2mat(t1) >= spectrum(1);
    
    t2 = {McNamaraBoswell(:).maxEx};
    t2indx = cellfun(@isempty,t2);
    t2(t2indx) = {0};
    maxEx = cell2mat(t2) >= spectrum(1);
    
    
    IDs = IDs(stokes & (maxAbs | maxEx));
    nExamples = length(IDs);
    indx = randperm(nExamples);
    
    % If we happened to request more fluorescent compounds than available,
    % we just provide what's available
    IDs = IDs(indx(1:min(nCompounds,nPatches)));
    
else
    IDs = compoundIDs;
end

nCompounds = length(IDs);
if isscalar(qe), qe = qe*ones(nCompounds,1); end

emission = zeros(length(McNamaraBoswell(1).wave),nCompounds);
absorption = zeros(length(McNamaraBoswell(1).wave),nCompounds);
names = cell(nCompounds,1);

for i=1:nCompounds
    indx = IDs(i);
    
    emission(:,i) = McNamaraBoswell(indx).emission;
    
    % If we don't have absorption data we use excitation
    if ~isempty(McNamaraBoswell(indx).absorption)
        absorption(:,i) = McNamaraBoswell(indx).absorption;
    else
        absorption(:,i) =  McNamaraBoswell(indx).excitation;
    end
    
    names{i} = [McNamaraBoswell(indx).dye ' ' McNamaraBoswell(indx).solvent];
end

% Normalize curves
emission = interp1(McNamaraBoswell(1).wave,emission,spectrum);
emission = emission*diag(1./max(emission));

absorption = interp1(McNamaraBoswell(1).wave,absorption,spectrum);
absorption = absorption*diag(1./(sum(absorption)*deltaL));

if isempty(weights)
    % Create the mixing matrix
    weights = rand(nCompounds,nPatches);
    weights(weights <= 0.5) = 0;
    if size(weights,1) >= 2
        weights = weights*diag(1./sum(weights));
        
        % Check to make sure that we don't have NaN's
        weights(:,isnan(sum(weights))) = 1;
        weights = weights*diag(1./sum(weights));
    end
end


photons = zeros(nWaves,nPatches);
for i=1:nPatches
    
    nAbsorbedPhotons = absorption'*illuminantGet(ill,'photons')*deltaL;
    emitted = emission*diag(nAbsorbedPhotons)*diag(qe)*weights(:,i);
    photons(:,i) = emitted;
end

photons = reshape(photons',h,w,nWaves);

photons = imageIncreaseImageRGBSize(photons,patchSize);

% Make it the right patch size
% fluorescenceChartObject.data.photons = photons;
fChart = fluorescenceSceneSet(fChart,'photons',photons);

end
