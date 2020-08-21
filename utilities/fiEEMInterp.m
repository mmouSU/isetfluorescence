function result = fiEEMInterp(eem, varargin)
% Interpolate the Excitation and Emission Matrix along a certain direction.
% For now I assume the x dimension is excitation and y dimension is
% emission.
%% Parse input
varargin = ieParamFormat(varargin);
p = inputParser;
p.addRequired('eem', @isnumeric);
p.addParameter('oldwave', 400:10:700, @isnumeric);
p.addParameter('newwave', 400:10:700, @isnumeric);
p.addParameter('dimension', 'excitation',@ischar);

p.parse(eem, varargin{:});
eem = p.Results.eem;
oldWave = p.Results.oldwave;
newWave = p.Results.newwave;
dimension = p.Results.dimension;

%% Do interpolation

switch dimension
    case 'excitation'
        tmp = interp1(oldWave, eem', newWave);
        result = tmp';
    case 'emission'
        result = interp1(oldWave, eem, newWave);
    case 'both'
        tmp = fiEEMInterp(eem, 'old wave',oldWave, 'new wave', newWave,...
                                'dimension', 'excitation');
        result = fiEEMInterp(tmp, 'old wave', oldWave, 'newwave', newWave,...
                                'dimension', 'emission');
    otherwise
        error('Unkown dimension %s', dimension);
end
result(isnan(result)) = 0;
if ~isequal(max(result(:)),1)
    warning('Peak excitation different from 1, rescaling')
    result = result/max(result(:));
end
%{
ieNewGraphWin;
imagesc(result)
%}
end