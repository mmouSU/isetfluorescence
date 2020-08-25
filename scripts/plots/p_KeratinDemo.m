%% Plot Keratin EEM

theseWaves = 300:5:705;
chdir(fullfile(fiToolboxRootPath,'data','Keratin'));
Keratin = fluorophoreRead('Keratin.mat','wave',theseWaves);
eemEnergy = Keratin.eemenergy;
fig = ieNewGraphWin;
clim = [0 1];
imagesc(theseWaves, theseWaves, eemEnergy, clim);
colorbar
xticks(300:50:700)
yticks(300:50:700)
xlabel('Excitation wavelength (nm)')
ylabel('Emission wavelength (nm)')
thisAxis = gca;
thisAxis.LineWidth = 2;
thisAxis.FontSize = 16;
set(gca,'GridLineStyle','--')

%% Plot several emission spetra under different wavelength
sampleWave = [375; 405; 435];
ieNewGraphWin;
hold all
for ii=1:numel(sampleWave)
    thisEmission = eemEnergy(:, theseWaves==sampleWave(ii));
    plot(theseWaves, thisEmission)
end
legend('375 nm', '405 nm', '435 nm');
title('Keratin Emission at different excitation wavelength');
ax = gca;
ax.FontSize=16;
xlabel('Wavelength (nm)');
ylabel('Intensity (a.u.)');