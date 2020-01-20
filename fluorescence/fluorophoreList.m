function fluorophores = fluorophoreList(fnames,wave)
% Return a cell array of fluorophores based on fnames
%
% Syntax
%   fluorophores = fluorophoreList(fnames,wave)
%
% Brief description
%
% Input
%  fnames:  Cell array of fluorophore names
%  wave:    Wavelength samples
%
% Optional key/value pairs
%
% Returns
%
%
% See also
%

%% Input parsing

%%
fluorophores = cell(numel(fnames),1);
for ii=1:numel(fnames)
    fluorophores{ii} = fiReadFluorophore(fnames{ii},'wave',wave);
end

