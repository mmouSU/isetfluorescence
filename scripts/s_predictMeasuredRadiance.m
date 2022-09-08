%% s_predictMeasuredRadiance
%
% Predict the radiance that we should measure of a slide with purported
% fluorophore when illuminated with a light
%
% To make this prediction, we use the following information
% 1. The excitation and emission spectra of the fluorophore
%       The excitation spectral curve (vector) describes how effective each wavelength is at
%       generating the emission
%       The emission spectral curve (vector) describes the spectral energy in the emitted light.
% 2. The spectral energy in the illuminating light
% 3. The spectral transmittance of the filter placed on the spectroradiometer
%       The fluorescence is small relative to the energy in the reflected
%       light, and the spetroradiometer cannot measure both the reflected light
%       and the fluorescence.  Hence, we use a Y44 filter to block light that is
%       great then 425 nm.  Presumably, this also blocks the Raman scattering of
%       the light as well
%
% To calculate the predicted radiance
% 1. Read in a fluorophore and create an excitation emission matrix (EEM)
% 2. Read in the light and take the product of the EEM and the light to get
%       the expected radiance
% 3. Multiply the expected radiance with the longpass filter to remove
% light that is not passed to the spectrophotometer
%
% J. Farrell December 2019
%
% See also
%    s_flSceneExample, s_CameraFluorescence, s_sceneFluorescence
%

%%
cd /users/joyce/Github/isetcam/
addpath(genpath(pwd))
cd ../isetfluorescence/
addpath(genpath(pwd))

ieInit

%% Select the range of wavelengths over which this calculation will be made
theseWaves = 300:5:800;
nWaves = numel(theseWaves);

%%  1. Read in a fluorophore and create an excitation emission matrix (EEM)
chdir(fullfile(fiToolboxRootPath,'data','webfluor'));

% Remove collagen and NADH and see the difference in the predictions (should be very little)
fluorophoreList = {'NADH_webfluor.mat','elastin_webfluor.mat','collagen1.mat','FAD_webfluor', 'protoporphyrin'};
dMatrix = zeros(nWaves,nWaves,numel(fluorophoreList));

for ii=1:numel(fluorophoreList)
    fName           = fullfile(fiToolboxRootPath,'data','webfluor',fluorophoreList{ii});
    fluor(ii)       = fiReadFluorophore(fName,'wave',theseWaves); 
    dMatrix(:,:,ii) = fluorophoreGet(fluor(ii),'eem');
end

%% 2. Read in the light 
illName = {'light405.mat','light450.mat' }; % repository has to be on the path
ill = zeros(nWaves,numel(illName));
for ii=1:length(illName)
    tmp = ieReadSpectra(illName{ii},theseWaves);
    tmp = mean(tmp,2);
    ill(:,ii) = tmp/max(tmp(:));
end
%{
 ieNewGraphWin
plot(theseWaves,ill);
 hold on; grid on; xlabel('Wavelength (nm)');ylabel('Relative radiance');
 legend(illName,'FontSize',18)
%}
%% and take the product of the EEM and the light to get the expected radiance

radiance = zeros(nWaves,length(illName),numel(fluorophoreList));
for ii =1:length(illName)
    for jj=1:numel(fluorophoreList)
        radiance(:,ii,jj) = dMatrix(:,:,jj)*ill(:,ii);
    end
end

%{
 for ii = 1:length(illName)
     ieNewGraphWin;
     for jj = 1:numel(fluorophoreList)
         plot(theseWaves,radiance(:,ii,jj));hold on;
     end
     grid on; xlabel('Wavelength (nm)');ylabel('Relative radiance');
     title(illName(ii));
     legend(fluorophoreList,'FontSize',18);
     ylim([0 4]);
end
%}
%% 3. Multiply the expected radiance with the longpass filter to remove light that is not passed to the spectrophotometer
% Remember that this prediction is assuming that there is no reflected light
% Place a shortpass filter on the light to minimize reflected light

% LPfilter = ieReadSpectra('Y44.mat',theseWaves); % longpass with 425 nm
% cutoff ... similar to the filter that is on the OralEye camera
LPfilter = ieReadSpectra('HoyaK2.mat',theseWaves); % Longpass with 475 nm cutoff
% we will use a 475 longpass filter
ieNewGraphWin;
plot(theseWaves,LPfilter);

for ii = 1:length(illName)
    ieNewGraphWin;
    for jj = 1:numel(fluorophoreList)
        plot(theseWaves, LPfilter .* radiance(:,ii,jj),'linewidth',2); hold on;
    end
    grid on; xlabel('Wavelength (nm)');ylabel('Relative radiance');
    title(illName(ii));
    legend(fluorophoreList,'FontSize',18);
    ax = gca;
    ax.FontSize=16;
    
end

% For each illuminant, add the radiance across fluorophores in order to get
% total predicted radiance assuming all fluorophores are stimulated
TotalRadiance = LPfilter.* sum(radiance(:,:,:),3);
for ii = 1:length(illName)
    ieNewGraphWin;
    plot(theseWaves,TotalRadiance(:,ii),'k','linewidth',3);
    title(illName(ii));
    ax = gca;
    ax.FontSize=16;
    ylim([0 4]);
end


%% without the LP filter
% For each illuminant, add the radiance across fluorophores in order to get
% total predicted radiance assuming all fluorophores are stimulated
TotalRadiance = sum(radiance(:,:,:),3);
for ii = 1:length(illName)
    ieNewGraphWin;
    plot(theseWaves,TotalRadiance(:,ii),'k','linewidth',3);
    title(illName(ii));
    ax = gca;
    ax.FontSize=16;
    ylim([0 4]);
end

%% END

