function [ scene ] = fiSceneAddFluorophore( scene, fl )

illSpd = sceneGet(scene,'illuminant photons');
wave = sceneGet(scene,'wave');
sz = sceneGet(scene,'size');
rePhotons = sceneGet(scene,'photons');

fl = fluorophoreSet(fl,'wave',wave);

ex = fluorophoreGet(fl,'excitation photons');
em = fluorophoreGet(fl,'emission photons');
qe = fluorophoreGet(fl,'qe');
deltaL = fluorophoreGet(fl,'deltaWave');

DM = tril(qe*deltaL*em*ex',-1);

flSpd = DM*illSpd;


flPhotons = repmat(shiftdim(flSpd,-2),[sz(:); 1]);

scene = sceneSet(scene,'photons',rePhotons + flPhotons);



end

