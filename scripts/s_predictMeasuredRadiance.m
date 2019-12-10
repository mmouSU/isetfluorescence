%% s_predictMeasuredRadiance

% The purpose of this script is to predict the radiance that we should
% measure of a slide with purported fluorophore when illuminated with a
% light

% To make this prediction, we use the following information
% 1. We need to know the excitation and emission spectra of the fluorophore
%       The excitation spectral curve (vector) describes how effective each wavelength is at
%       generating the emission
%       The emission spectral curve (vector) describes the spectral energy in the emitted light.
% 2. We need to know the spectral energy in the illuminating light
% 3. We need to know the spectral transmittance of the filter placed on the
% spectroradiometer
%       The fluorescence is small relative to the energy in the reflected
%       light, and the spetroradiometer cannot measure both the reflected light
%       and the fluorescence.  Hence, we use a Y44 filter to block light that is
%       great then 425 nm.  Presumably, this also blocks the Raman scattering of
%       the light as well

% To calculate the predicted radiance
% 1. Read in a fluorophore and create an excitation emission matrix (EEM)
% 2. Read in the light and take the product of the EEM and the light to get
%       the expected radiance
% 3. Multiply the expected radiance with the longpass filter to remove
% light that is not passed to the spectrophotometer

% J. Farrell December 2019

%% Select the range of wavelengths over which this calculation will be made
theseWaves = 300:5:800;
nWaves = numel(theseWaves);

%%  1. Read in a fluorophore and create an excitation emission matrix (EEM)
fluorophoreList = {'NADH.mat','elastin.mat','collagen1.mat'};
dMatrix = zeros(nWaves,nWaves,numel(fluorophoreList));

for ii=1:numel(fluorophoreList)
    fName           = fullfile(fiToolboxRootPath,'data','webfluor',fluorophoreList{ii});
    fluor(ii)       = fiReadFluorophore(fName,'wave',theseWaves); 
    dMatrix(:,:,ii) = fluorophoreGet(fluor(ii),'donaldson matrix');
end

%% 2. Read in the light and take the product of the EEM and the light to get the expected radiance
illName = {'OralEye_385.mat'}; % oraleye repository has to be on the path
ill = zeros(nWaves,numel(illName));
for ii=1:length(illName)
    tmp = ieReadSpectra(illName{ii},theseWaves);
    tmp = mean(tmp,2);
    ill(:,ii) = tmp/max(tmp(:));
end
% 
% ieNewGraphWin
% plot(theseWaves,ill);
% hold on; grid on; xlabel('Wavelength (nm)');ylabel('Relative radiance');
% legend(illName,'FontSize',18)

radiance = zeros(nWaves,numel(fluorophoreList));
for jj=1:numel(fluorophoreList)
    radiance(:,jj) = dMatrix(:,:,jj)*ill(:);
end

ieNewGraphWin
plot(theseWaves,radiance);
hold on; grid on; xlabel('Wavelength (nm)');ylabel('Relative radiance');
legend(fluorophoreList,'FontSize',18)


%% 3. Multiply the expected radiance with the longpass filter to remove light that is not passed to the spectrophotometer

LPfilter = ieReadSpectra('Y44.mat',theseWaves);
filteredRadiance = LPfilter .* radiance;

ieNewGraphWin
plot(theseWaves,filteredRadiance);
hold on; grid on; xlabel('Wavelength (nm)');ylabel('Filtered radiance');
legend(fluorophoreList,'FontSize',18)
