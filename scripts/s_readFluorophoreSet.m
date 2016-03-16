ieInit;

dirPath = fullfile(fiToolboxRootPath,'data','LifeTechnologies');

[flSet, ids] = fiReadFluorophoreSet(dirPath,'peakExRange',[420 Inf],...
                                            'peakEmRange',[0 680],...
                                            'stokesShiftRange',[20 40]);