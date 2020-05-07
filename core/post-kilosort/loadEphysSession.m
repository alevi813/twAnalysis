function [PDS, flagBits, flagData, flagPlexonTime] = loadEphysSession(pdsfile, plxname)

% INPUTS
% pdsfile (path and name)
% plxname (path and name)

% pdsfile = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160615twag.twag1401.PDS';
% plxname = '/Users/Aaron/Dropbox/nancy_twag_20160615_exp1.plx';

% pdsfile = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160622twag.twag1345.PDS';
% plxname = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/ephys/nTwag20160622exp1.plx';

% pdsfile = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160726twag.twag1540.PDS';
% plxname = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/ephys/nTwag20160726_mt_l5p4_exp1.plx';
% ch 10

% pdsfile = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160729twag.twag1317.PDS';
% plxname = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/ephys/nTwag_20160729_mt_l6p4_exp2.plx';
% % 1-3, 9-22 ---- something's wrong with CP here...

pdsfile = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160802twag.twag1356.PDS';
plxname = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/ephys/nTwag_20160802_l5p3_exp1.plx';
% 5 7-16

% pdsfile = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160803twag.twag1233.PDS';
% plxname = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/ephys/nTwag_20160803_l4p4_exp1.plx';
% 11 - 16 

% pdsfile = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160804twag.twag1441.PDS';
% plxname = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/ephys/nTwag_20160804_l6p4_exp1.plx';
%9-12 good, 11 really really good

%%% FREE CHOICE
% pdsfile = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160803twag.FreeChoice1528.PDS';
% plxname = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/ephys/nFreeChoice_20160803_l4p4_exp1.plx';

load(pdsfile, '-mat')
plx=ephys.readPlx(plxname, true);

b=plx.eventChannels.values;
b=double(typecast(int16(b),'uint16'));

flagBits = log(bitshift(b,-8))/log(2) + 1;
flagData = mod(b,2^8);

flagPlexonTime = double(plx.eventChannels.events);
%[PL2PTBfit,PL2PTB,PTB2PL, maxreconstructionerror, flagPlexonTrialStartTime ] = ephys.syncPlexonClock(PDS, plxname);
end