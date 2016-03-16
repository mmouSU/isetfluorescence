function [fChart, excitation, emission, names, fluorophoreIDs, donaldsonMat] = fluorescenceChartCreate(patchSize,h,w,spectrum,ill,varargin)

p = inputParser;
p.KeepUnmatched = true;
p.addParamValue('qe',1,@isfloat);
p.addParamValue('nFlPerPatch',1,@isfloat);
p.addParamValue('minStokes',0,@isfloat);
p.addParamValue('maxStokes',Inf,@isfloat);
p.addParamValue('fluorophoreIDs',{});

p.parse(varargin{:});
inputs = p.Results;

% This is the size in pixels of each testchart patch
if ieNotDefined('patchSize'), patchSize = 16;   end

% Define the size of the chart
if ieNotDefined('h'), h = 4; end
if ieNotDefined('w'), w = 6; end


deltaL = spectrum(2) - spectrum(1);
nWaves = length(spectrum);
nPatches = h*w;
fName = fullfile(isetRootPath,'data','surfaces','McNamaraBoswell.mat');
flDataset = load(fName);

excitation = cell(h,w);
emission = cell(h,w);
names = cell(h,w);
fluorophoreIDs = cell(h,w);
donaldsonMat = cell(h,w);
flRadiance = zeros(h,w,nWaves);


fChart = fluorescenceSceneCreate('Fluorescence chart',spectrum);


% If we don't specify patches, we draw them at random from the databse, so
% that they satisfy the conditions.
if isempty(inputs.fluorophoreIDs)
    
    % Select all the fluorophores that satisfy the conditions
    
    % 1. They fall within the Stokes limits
    IDs = 1:length(flDataset.McNamaraBoswell);
    stokes = ([flDataset.McNamaraBoswell(:).stokesShift] >= inputs.minStokes) & ...
              ([flDataset.McNamaraBoswell(:).stokesShift] < inputs.maxStokes);
    
    % 2. Restrict the absorption and emission peak to be within the
    % discretized range of values
    t1 = {flDataset.McNamaraBoswell(:).maxAbs};
    t1indx = cellfun(@isempty,t1);
    t1(t1indx) = {0};
    maxAbs = cell2mat(t1) >= spectrum(1);
    
    t2 = {flDataset.McNamaraBoswell(:).maxEx};
    t2indx = cellfun(@isempty,t2);
    t2(t2indx) = {0};
    maxEx = cell2mat(t2) >= spectrum(1);
    
    
    IDs = IDs(stokes & (maxAbs | maxEx));
    nExamples = length(IDs);
    indx = randperm(nExamples);
    
    % If we still have more examples in the database we sample without
    % repetition
    while length(indx) < nPatches*inputs.nFlPerPatch
        indx = [indx, indx];
    end
    
    for ww=1:w
        for hh=1:h
            ll = inputs.nFlPerPatch*h*(ww-1) + inputs.nFlPerPatch*(hh-1) + 1;
            ul = inputs.nFlPerPatch*h*(ww-1) + inputs.nFlPerPatch*(hh);
            fluorophoreIDs{hh,ww} = IDs(indx(ll:ul));
        end
    end
    
else
    fluorophoreIDs = inputs.fluorophoreIDs;
end


for ww=1:w
   for hh=1:h
      
       em = [];
       ex = [];
       nm = {};
       for i=1:length(fluorophoreIDs{hh,ww})
           
           indx = fluorophoreIDs{hh,ww}(i);
           
           em = [em, flDataset.McNamaraBoswell(indx).emission];
           if ~isempty(flDataset.McNamaraBoswell(indx).absorption)
               ex = [ex, flDataset.McNamaraBoswell(indx).absorption];
           else
               ex = [ex, flDataset.McNamaraBoswell(indx).excitation];
           end
           
           nm{i} = [flDataset.McNamaraBoswell(indx).dye ' ' flDataset.McNamaraBoswell(indx).solvent];
           
       end
       
       em = interp1(flDataset.McNamaraBoswell(1).wave,em,spectrum);
       ex = interp1(flDataset.McNamaraBoswell(1).wave,ex,spectrum);
       
       ex = ex*diag(1./max(ex));
       em = em*diag(1./sum(em*deltaL));
       
       % nAbsorbedPhotons = inputs.qe*ex'*illuminantGet(ill,'photons')*deltaL;
     
       emission{hh,ww} = em;
       excitation{hh,ww} = ex;
       names{hh,ww} = nm;
       
       % The Donaldson matrix is normalized so that donaldsonMat*illuminant
       % = fluorescentRadiance
       fullDonaldson = inputs.qe*em*ex'/length(fluorophoreIDs{hh,ww});
       donaldsonMat{hh,ww} = tril(fullDonaldson,-1);
       
       % figure; imagesc([fullDonaldson donaldsonMat{hh,ww}]);
       
       % flRadiance(hh,ww,:) = donaldsonMat{hh,ww}*illuminantGet(ill,'photons')*deltaL;
       
   end
end

% photons = imageIncreaseImageRGBSize(flRadiance,patchSize);

% Make it the right patch size
% fluorescenceChartObject.data.photons = photons;
% fChart = fluorescenceSceneSet(fChart,'photons',photons);

fChart = fluorescenceSceneFromCell(emission,excitation,patchSize,spectrum,inputs.qe,ill);

end
