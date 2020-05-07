%%%% preanalysis pipeline

% save a pp file
baseDir = '/Volumes/HukLab/Macaque/Projects/twag/leo/';
pdsf    = '20191218/leo20191218_l3p2_final/leo20191218twag.twag1255.PDS';

ppNum   = 35;

pp = analyzeBehaviorSession(baseDir, pdsf, 1, 1, ppNum);


%% twagTrimmedDataFactory

load([baseDir pdsf], '-mat');

plxfile = '20191218/leo20191218_l3p2_final/leo20191218_l3p2.pl2';
plxfile = [baseDir plxfile];

syncClock = ephys.syncPlexonClock(PDS, plxfile);


save([baseDir '20191218/leo20191218_l3p2_final/syncClock.mat'], 'syncClock');

%% twagEphysAnalysis ---- bp after twagPSTH

%% twagEphysAnalysis_postKilo --- +phy

%% twagEphysAnalysis_postKilo -- run with bp for goodClust thresholding

%% run pds2stim and kilo2n