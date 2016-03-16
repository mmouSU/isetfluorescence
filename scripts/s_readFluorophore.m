ieInit;

fName = fullfile(fiToolboxRootPath,'data','LifeTechnologies','AlexaFluor555');

fl = fiReadFluorophore(fName);

wave = fluorophoreGet(fl,'wave');
figure; plot(wave,[fluorophoreGet(fl,'excitation') fluorophoreGet(fl,'emission')]);

fl = fiReadFluorophore(fName,'wave',400:10:700);
wave = fluorophoreGet(fl,'wave');
figure; plot(wave,[fluorophoreGet(fl,'norm excitation') fluorophoreGet(fl,'norm emission')]);

