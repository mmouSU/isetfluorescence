function fChart = fluorescenceSceneFromCell(emissionCell,excitationCell,patchSize,wave,qe,illuminant)

fChart = fluorescenceSceneCreate('Fluorescence chart',wave);

% fluorescenceChartObject.name = 'Fluorescence chart';
% fluorescenceChartObject.type = 'fluorescence scene';

% This is the size in pixels of each testchart patch
if ieNotDefined('patchSize'), patchSize = 16;   end

% This parameter specifies the quantuum efficiency, i.e. the ratio between
% the number of emited to absorbed photons.
if ieNotDefined('qe'), qe = 1; end

if ieNotDefined('illuminant')
    illuminant = illuminantCreate('equalEnergy',[],100,wave)';
end
fChart = fluorescenceSceneSet(fChart,'illuminant',illuminant);

fChart = fluorescenceSceneSet(fChart,'excitation',excitationCell);
fChart = fluorescenceSceneSet(fChart,'emission',emissionCell);
fChart = fluorescenceSceneSet(fChart,'patchSize',patchSize);
fChart = fluorescenceSceneSet(fChart,'quantumEfficiency',qe);

end
