% This script demonstrates the use of functions associated with creating, and 
% operating on a fluorescent scene in the fiToolbox.
%
% Copyright, Henryk Blasinski 2016.

close all;
clear variables;
clc;

ieInit;

wave = 400:10:700;

% First let's create a standard reflective ISET scene, it will be used as the
% reflective component in subsequent tests. 
scene = sceneCreate('macbeth d65',wave);
ieAddObject(scene);
sceneWindow();

% Create a scene with uniform fluorescent properties of the 50th compound
% from the McNamara-Boswell data set. Set the practical quantum efficiency
% to 0.5
singleFlScene = fluorescentSceneCreate('type','onefluorophore','fluorophoreIDs',50);
singleFlScene = fluorescentSceneSet(singleFlScene,'qe',0.5);

% Combine the reflective and fluorescent scenes.
sceneWithSingleFl = fiSceneAddFluorescence(scene,singleFlScene);
sceneWithSingleFl = sceneSet(sceneWithSingleFl,'name','Scene with one fluorophore, D65');
ieAddObject(sceneWithSingleFl);
sceneWindow();

% To investigate fluorescence under different illuminant we first need to
% change the ISET scene illuminant, and then add the fluorescent scene.
ill = illuminantCreate('blackbody',wave,2000);
scene2 = sceneAdjustIlluminant(scene,illuminantGet(ill,'energy'));
scene2 = sceneSet(scene2,'name','Macbeth, 2000K');
ieAddObject(scene2);

scene2WithSingleFl = fiSceneAddFluorescence(scene2,singleFlScene);
sceneWithSingleFl = sceneSet(sceneWithSingleFl,'name','Scene with one fluorophore, 2000K');
ieAddObject(sceneWithSingleFl);
sceneWindow();


% Let's create a 2 x 3 multi-fluorophore scene with two fluorophores per
% spatial location and some restrictions on the choice of fluorophores.
% When we combine this fluorescent scene with a Macbeth chart, fluorescence
% properties will be shared between four reflectance patches.

multiFlScene = fluorescentSceneCreate('height',2,'width',3,'nFluorophores',1,...
                                      'stokesShiftRange',[20 Inf],...
                                      'peakEmRange',[450 550],...
                                      'wave',wave);
                                  
sceneWithMultiFl = fiSceneAddFluorescence(scene,multiFlScene);
sceneWithMultiFl = sceneSet(sceneWithMultiFl,'name','Scene with fluorophore array');
ieAddObject(sceneWithMultiFl);
sceneWindow();

% We can extract the reference fluorescence spectra and plot them for every
% patch. Note that fluorescence properties are shared between four patches.
emRef = fluorescentSceneGet(multiFlScene,'emission reference','sceneSize',[4 6]);
exRef = fluorescentSceneGet(multiFlScene,'excitation reference','sceneSize',[4 6]);

figure;
for x=1:6
    for y=1:4
    
        patchID = (x-1)*4+y;
        figID = (y-1)*6 + x;
        
        subplot(4,6,figID);
        hold on; grid on; box on;
        
        fl = [emRef(:,patchID) exRef(:,patchID)];
        
        plot(wave,fl*diag(1./max(fl)),'lineWidth',1.5);
        xlim([min(wave) max(wave)]);
        ylim([-0.05 1.05]);
        xlabel('Wavelength, nm');
    end
end


