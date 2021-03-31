% Example of how to get EEMs and visualize them

%%
ieInit;
%%
wave = 365:5:705; % Specify a wavelength range
%%
flNameList = {'FAD_webfluor', 'NADH_webfluor', 'collagen1', 'elastin_webfluor',...
               'Porphyrins'}; % Fluorophores used from database
eemList = cell(1, numel(flNameList));

for ii=1:numel(flNameList)
    fluorophore = fluorophoreRead(flNameList{ii}, 'wave', wave); % Returns a fluorophore struct
    eemList{ii} = fluorophoreGet(fluorophore, 'eem energy normalize'); % Generate normalized eem
end

%% Visualize eems

for ii=1:numel(flNameList)
    ieNewGraphWin;
    imagesc(wave, wave, eemList{ii}');
    title(flNameList{ii});
    
    xticks(350:50:700)
    yticks(350:50:700)
    ylabel('Excitation wavelength (nm)')
    xlabel('Emission wavelength (nm)')
    colorbar
    
    thisAxis = gca;
    thisAxis.LineWidth = 2;
    set(gca,'GridLineStyle','--')
    set(gca,'FontSize',20)
    set(gca, 'FontWeight', 'bold')
    set(gca,'YDir','normal')
    axis square
end