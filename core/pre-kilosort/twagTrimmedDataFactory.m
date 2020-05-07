%function [pdsfile, plxfile, PDS, ppNum, PL2PTBfit, PL2PTB, PTB2PL, maxreconstructionerror, flagPlexonTrialStartTime] = twagTrimmedDataFactory(subject, eList)
function [pdsfile, plxfile, PDS, ppNum, syncClock] = twagTrimmedDataFactory(subject, eList)

if strcmp(subject, 'leo')
    %     dirDataPDS = '/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/leo/';
    %     dirDataPlx = '/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/leo/';
    dirDataPDS = '/Volumes/HukLab/Macaque/Projects/twag/leo/';
    dirDataPlx = '/Volumes/HukLab/Macaque/Projects/twag/leo/';
else
%     dirDataPDS   = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/';
%     dirDataPlx   = '/Volumes/LKCLAB/EPHYS/DATA/PLXdata/';
    dirDataPDS = '/Volumes/HukLab/Macaque/Projects/twag/nancy/';
    dirDataPlx = '/Volumes/HukLab/Macaque/Projects/twag/nancy/';
end

nExp         = length(eList);
ppAll        = cell(nExp,1);

for iE = 1:nExp
    expID = eList(iE);
    
    switch subject
        
        case 'nancy'
            switch expID
                
                case 1 % ch 10 Lpref ++++++++++++++++++++++++++++++++++++++
                    pdsfile = 'nancy/nancy20160726twag.twag1540.PDS';
                    plxfile = 'nTwag20160726_mt_l5p4_exp1.plx';
                    ppNum   = 58;
                    
                case 2 % 5, 7-16 ++++++++++++++++++++++++++++++++++++++++++
                    pdsfile = 'nancy/nancy20160802twag.twag1356.PDS';
                    plxfile = 'nTwag_20160802_l5p3_exp1.plx';
                    ppNum   = 60;
                case 3 % 11 - 16 ++++++++++++++++++++++++++++++++++++++++++
                    pdsfile = 'nancy/nancy20160803twag.twag1233.PDS';
                    plxfile = 'nTwag_20160803_l4p4_exp1.plx';
                    ppNum   = 61;
                case 4 %9 - 12  +++++++++++++++++++++++++++++++++++++++++++
                    pdsfile = 'nancy/nancy20160804twag.twag1441.PDS';
                    plxfile = 'nTwag_20160804_l6p4_exp1.plx';
                    ppNum   = 62;
                    
                case 5 %       ++++++++++++++++++++++++++++++++++++++++++++
                    %3 - 7, 8-12 are there but not great, 13-17 are there but weird and Lpref
                    pdsfile = 'nancy/nancy20160812twag.twag1151.PDS';
                    plxfile = 'nTwag_2016081216_exp1.plx';
                    ppNum   = 64;
                case 6 % 7,8 ok Lpref; 13-18 seem good, but DS switches from L to R at 17?
                    % 19 & 20 ok? this is a questionable data set.
                    pdsfile = 'nancy/nancy20160823twag.twag1405.PDS';
                    plxfile = 'nTwag_20160823_l4p4_exp1.plx';
                    ppNum = 65;
                case 7 % ++++++++++++ 1-24 ok (1,2 best DS) 11-13 weirdly most response to dist 3
                    % 14 a little weird looking, 5-19 great, 23-24 great
                    pdsfile = 'nancy/nancy20160825twag.twag1324.PDS';
                    plxfile = 'nTwag_20160825_l4p4_exp2.plx';
                    ppNum   = 66;
                case 8 % ++++++++++++++++++++++++ 1-6 (3,5,6 best) 7-18
                    pdsfile = 'nancy/nancy20160831twag.twag1330.PDS';
                    plxfile = 'nTwag_20160831_l4p4_exp1.plx';
                    ppNum   = 67;
                case 9 % ++++++++++++ 3-12  (3, 4, 5, 6, 12 best) 20, 21 ok. Some other weird ones that could be alright.
                    pdsfile = 'nancy/nancy20160906twag.twag1421.PDS';
                    plxfile = 'nTwag20160905_l4p4_exp1.plx';
                    ppNum   = 68;
                case 10 % ++++++++++++ 1-5 good
                    pdsfile = 'nancy/nancy20160908twag.twag1248.PDS';
                    plxfile = 'nTwag_20160908_l4p4_exp1.plx';
                    ppNum   = 69;
                case 11 %%% LATE STIMULUS %%% ++++++++
                    pdsfile = 'nancy/nancy20160909twag.twag1247.PDS';
                    plxfile = 'nTwag_20160909_l4p4_exp1.plx';
                    ppNum   = 70;
                case 12 % ++++++++
                    pdsfile = 'nancy/nancy20160913twag.twag1452.PDS';
                    plxfile = 'nTwag_20160913_l4p4_exp1.plx';
                    ppNum   = 71;
                case 13 % ++++++++
                    pdsfile = '20160915/nancy20160915twag.twag1328.PDS';
                    plxfile = '20160915/nTwag20160915_l4p4_exp1.plx';
                    ppNum   = 72;
                case 14 % ++++++++
                    pdsfile = 'nancy/nancy20160920twag.twag1310.PDS';
                    plxfile = 'nTwag20160920_l4p4_exp1.plx';
                    ppNum   = 73;
                case 15 % ++++++++ ch 19, 20 good. that's about it
                    pdsfile = 'nancy/nancy20160921twag.twag1449.PDS';
                    plxfile = 'nTwag_20160921_l4p4_exp1.plx';
                    ppNum   = 74;
                case 16 % ++++++++
                    pdsfile = 'nancy/nancy20160923twag.twag1314.PDS';
                    plxfile = 'nTwag20160923_l4p4_exp1.plx';
                    ppNum   = 75;
                case 17 % ++++++++
                    pdsfile = 'nancy/nancy20160928twag.twag1401.PDS';
                    plxfile = 'nTwag_20160928_l4p4_exp1.plx';
                    ppNum   = 76;
                case 18 % ++++++++
                    pdsfile = 'nancy/nancy20161005twag.twag1454.PDS';
                    plxfile = 'nTwag_20161005_l5p3_exp1.plx';
                    ppNum   = 77;
                case 19 % ++++++++ 21-24 are great... rest are ok/kinda weird
                    pdsfile = 'nancy/nancy20161013twag.twag1255.PDS';
                    plxfile = 'nTwag_20162013_l5p3_exp1.plx';
                    ppNum   = 78;
                case 20 % ++++++
                    pdsfile = 'nancy/nancy20161014twag.twag1359.PDS';
                    plxfile = 'nTwag_20161014_l5p3_exp2.plx';
                    ppNum   = 79;
                case 21 % +++++
                    pdsfile = 'nancy/nancy20161018twag.twag1330.PDS';
                    plxfile = 'nTwag_20162018_l5p3_exp1.plx';
                    ppNum   = 80;
                case 22 % start early +++++
                    pdsfile = 'nancy/nancy20161019twag.twag1403.PDS';
                    plxfile = 'nTwag_20161019_l5p3_exp1.plx';
                    ppNum   = 81;
                case 23 % +++++
                    pdsfile = 'nancy/nancy20161020twag.twag1301.PDS';
                    plxfile = 'nTwag_20161020_l5p3_exp1.plx';
                    ppNum   = 82;
                case 24 % ++++++
                    pdsfile = 'nancy/nancy20161021twag.twag1349.PDS';
                    plxfile = 'nTwag_20161021_l5p3_exp1.plx';
                    ppNum   = 83;
                    
                case 25 % ++++++
                    pdsfile = 'nancy/nancy20161027twag.twag1343.PDS';
                    plxfile = 'nTwag_20161027_l5p4_exp1.plx';
                    ppNum   = 85;
                case 26 % ++++++
                    pdsfile = 'nancy/nancy20161028twag.twag1315.PDS';
                    plxfile = 'nTwag_20161028_l5p4_exp1.plx';
                    ppNum   = 86;
                case 27 % ++++++
                    pdsfile = 'nancy/nancy20161101twag.twag1329.PDS';
                    plxfile = 'nTwag_20161101_l5p4_exp1.plx';
                    ppNum   = 87;
                case 28 % ++++++
                    pdsfile = 'nancy/nancy20161102twag.twag1230.PDS';
                    plxfile = 'nTwag_20161102_l5p4_exp1.plx';
                    ppNum   = 88;
                case 29 % ++++++
                    pdsfile = 'nancy/nancy20161103twag.twag1259.PDS';
                    plxfile = 'nTwag_20161103_l5p3_exp1.plx';
                    ppNum   = 89;
                case 30 % longer PDS
                    pdsfile = 'nancy/nancy20161109twag.twag1315.PDS';
                    plxfile = 'nTwag_20161109_l5p4_exp1.plx';
                    ppNum   = 90;
                case 31 % ++++ early cp session
                    pdsfile = 'nancy/nancy20161110twag.twag1339.PDS';
                    plxfile = 'nTwag_20161110_l5p4_exp1.plx';
                    ppNum   = 91;
                case 32 % ++++
                    pdsfile = 'nancy/nancy20161111twag.twag1257.PDS';
                    plxfile = 'nTwag_20161111_l5p4_exp1.plx';
                    ppNum   = 92;
                case 33 % ++++
                    pdsfile = 'nancy/nancy20161115twag.twag1427.PDS';
                    plxfile = 'nTwag_20161115_l5p4_exp1.plx';
                    ppNum   = 93;
                case 34 % ++++
                    pdsfile = 'nancy/nancy20161130twag.twag1345.PDS';
                    plxfile = 'nTwag_20161130_l5p4_exp1.plx';
                    ppNum   = 99;
                case 35 % ++++
                    pdsfile = 'nancy/nancy20161201twag.twag1407.PDS';
                    plxfile = 'nTwag_20161201_l5p4_exp1.plx';
                    ppNum   = 100;
                case 36 % ++++
                    pdsfile = 'nancy/nancy20161202twag.twag1433.PDS';
                    plxfile = 'nTwag_20161202_l5p4_exp1.plx';
                    ppNum   = 101;
                case 37 % ++++
                    pdsfile = 'nancy/nancy20161209twag.twag1246.PDS';
                    plxfile = 'nTwag_20161209_l5p4_exp1.plx';
                    ppNum   = 103;
            end
        case 'leo'
            switch expID
                case 1
                    pdsfile = '20190806/leo20190806_l3p2_2_final/leo20190806twag.twag1234.PDS';
                    plxfile = '20190806/leo20190806_l3p2_2_final/leo20190806_l3p2_2.pl2';
                    ppNum = 1;
                    
                case 2
                    pdsfile = '20190807/leo20190807_l4p2_2/leo20190807twag.twag1310.PDS';
                    plxfile = '20190807/leo20190807_l4p2_2/leo20190807_l4p2_2.pl2';
                    ppNum = 2;
                    
                case 3
                    pdsfile = '20190814/leo20190814_l4a1_1/leo20190814twag.twag1206.PDS';
                    plxfile = '20190814/leo20190814_l4a1_1/leo20190814_l4a1_1.pl2';
                    ppNum = 3;
                    
                case 4
                    pdsfile = '20190815/leo20190815_l3a1_1/leo20190815twag.twag1314.PDS';
                    plxfile = '20190815/leo20190815_l3a1_1/leo20190815_l3a1_1.pl2';
                    ppNum = 4;
                    
                case 5
                    pdsfile = '20190816/leo20190816_l4p2_1/leo20190816twag.twag1228.PDS';
                    plxfile = '20190816/leo20190816_l4p2_1/leo20190816_l4p2_1.pl2';
                    ppNum = 5;
                    
                case 6
                    pdsfile = '20190827/leo20190827_l3p1_final/leo20190827twag.twag1159.PDS';
                    plxfile = '20190827/leo20190827_l3p1_final/leo20190827_l3p1.pl2';
                    ppNum = 6;
                    
                case 7
                    pdsfile = '20190829/leo20190829_l4p2_1_final/leo20190829twag.twag1200.PDS';
                    plxfile = '20190829/leo20190829_l4p2_1_final/leo20190829_l4p2_1.pl2';
                    ppNum = 7;
                    
                case 8
                    pdsfile = '20190830/leo20190830_l5p2_final/leo20190830twag.twag1219.PDS';
                    plxfile = '20190830/leo20190830_l5p2_final/leo20190830_l5p2.pl2';
                    ppNum = 8;
                    
                case 9
                    pdsfile = '20190904/leo20190904_l4p1_final/leo20190904twag.twag1147.PDS';
                    plxfile = '20190904/leo20190904_l4p1_final/leo20190904_l4p1.pl2';
                    ppNum = 9;
                    
                case 10
                    pdsfile = '20190905/leo20190905_l3p2_final/leo20190905twag.twag1253.PDS';
                    plxfile = '20190905/leo20190905_l3p2_final/leo20190905_l3p2.pl2';
                    ppNum = 10;
                    
                case 11
                    pdsfile = '20190916/leo20190916_l4a1_final/leo20190916twag.twag1218.PDS';
                    plxfile = '20190916/leo20190916_l4a1_final/leo20190916_l4a1.pl2';
                    ppNum = 11;
                    
                case 12
                    pdsfile = '20190917/leo20190917_l4p1_final/leo20190917twag.twag1223.PDS';
                    plxfile = '20190917/leo20190917_l4p1_final/leo20190917_l4p1.pl2';
                    ppNum = 12;
                    
                case 13
                    pdsfile = '20190918/leo20910918_l3p1_final/leo20190919twag.twag1209.PDS';
                    plxfile = '20190918/leo20910918_l3p1_final/leo20910918_l3p1.pl2';
                    ppNum = 13;
                    
                case 14
                    pdsfile = '20191008/leo20191008_l4p2/leo20191008twag.twag1148.PDS';
                    plxfile = '20191008/leo20191008_l4p2/leo20191008_l4p2.pl2';
                    ppNum = 14;
                    
                case 15
                    pdsfile = '20191009/leo20191009twag.twag1203.PDS';
                    plxfile = '20191009/leo20191009_l4p1.pl2';
                    ppNum = 15;
                    
                case 16
                    pdsfile = '20191010/leo20191010_l3p1/leo20191010twag.twag1149.PDS';
                    plxfile = '20191010/leo20191010_l3p1/leo20191010_l3p1.pl2';
                    ppNum = 16;
                    
                case 17
                    pdsfile = '20191016/leo20191016_l4p1_final/leo20191016twag.twag1235.PDS';
                    plxfile = '20191016/leo20191016_l4p1_final/leo20191016_l4p1.pl2';
                    ppNum = 17;
                    
                case 18
                    pdsfile = '20191021/leo20191021_l4p1/leo20191021twag.twag1310.PDS';
                    plxfile = '20191021/leo20191021_l4p1/leo20191021_l4p1.pl2';
                    ppNum = 18;
                    
                case 19
                    pdsfile = '20191023/leo20191023_l4p1_final/leo20191023twag.twag1230.PDS';
                    plxfile = '20191023/leo20191023_l4p1_final/leo20191023_l4p1.pl2';
                    ppNum = 19;
                    
                case 20
                    pdsfile = '20191025/leo20191025_l3p1_final/leo20191025twag.twag1323.PDS';
                    plxfile = '20191025/leo20191025_l3p1_final/leo20191025_l3p1.pl2';
                    ppNum = 20;
                    
                case 21
                    pdsfile = '20191120/leo20191120_l4p1_final/leo20191120twag.twag1225.PDS';
                    plxfile = '20191120/leo20191120_l4p1_final/leo20191120_l4p1.pl2';
                    ppNum = 21;
                    
                case 22
                    pdsfile = '20191122/leo20191122_l4p2_final/leo20191122twag.twag1321.PDS';
                    plxfile = '20191122/leo20191122_l4p2_final/leo20191122_l4p2.pl2';
                    ppNum = 22;
                    
                case 23
                    pdsfile = '20191126/leo20191126_l4p1_final/leo20191126twag.twag1305.PDS';
                    plxfile = '20191126/leo20191126_l4p1_final/leo20191126_l4p1.pl2';
                    ppNum = 23;
                    
                case 24
                    pdsfile = '20191129/leo20191129_l4p2_final/leo20191129twag.twag1348.PDS';
                    plxfile = '20191129/leo20191129_l4p2_final/leo20191129_l4p2.pl2';
                    ppNum = 24;
                    
                case 25
                    pdsfile = '20191205/leo20191205_l3p2_final/leo20191205twag.twag1338.PDS';
                    plxfile = '20191205/leo20191205_l3p2_final/leo20191205_l3p2.pl2';
                    ppNum = 25;
                    
                case 26
                    pdsfile = '20191206/leo20191206_l5p2_final/leo20191206twag.twag1310.PDS';
                    plxfile = '20191206/leo20191206_l5p2_final/leo20191206_l5p2.pl2';
                    ppNum = 26;
                    
                case 27
                    pdsfile = '20191207/leo20191207_l3p2_final/leo20191207twag.twag1237.PDS';
                    plxfile = '20191207/leo20191207_l3p2_final/leo20191207_l3p2.pl2';
                    ppNum = 27;
                    
                case 28
                    pdsfile = '20191209/leo20191209_l5p2_final/leo20191209twag.twag1233.PDS';
                    plxfile = '20191209/leo20191209_l5p2_final/leo20191209_l5p2.pl2';
                    ppNum = 28;
                    
                case 29
                    pdsfile = '20191210/leo20191210_l3p2_final/leo20191210twag.twag1400.PDS';
                    plxfile = '20191210/leo20191210_l3p2_final/leo20191210_l3p2.pl2';
                    ppNum = 29;
                    
                case 30
                    pdsfile = '20191211/leo20191211_l3p2_final/leo20191211twag.twag1320.PDS';
                    plxfile = '20191211/leo20191211_l3p2_final/leo20191211_l3p2.pl2';
                    ppNum = 30;      
                    
                case 31
                    pdsfile = '20191212/leo20191212_l3p2_final/leo20191212twag.twag1243.PDS';
                    plxfile = '20191212/leo20191212_l3p2_final/leo20191212_l3p2.pl2';
                    ppNum = 31; 
                    
                case 32
                    pdsfile = '20191213/leo20191213_l5p2_final/leo20191213twag.twag1258.PDS';
                    plxfile = '20191213/leo20191213_l5p2_final/leo20191213_l5p2.pl2';
                    ppNum = 32;      
                                        
                case 33
                    pdsfile = '20191216/leo20191216_l3p2_final/leo20191216twag.twag1309.PDS';
                    plxfile = '20191216/leo20191216_l3p2_final/leo20191216_l3p2.pl2';
                    ppNum = 33;       
                    
                case 34
                    pdsfile = '20191217/leo20191217_l3p2_final/leo20191217twag.twag1432.PDS'; %%%% TWO PDS FILES (EYELINK CRASH)
                    plxfile = '20191217/leo20191217_l3p2_final/leo20191217_l3p2.pl2';
                    ppNum = 34;         
                    
                case 35
                    pdsfile = '20191218/leo20191218_l3p2_final/leo20191218twag.twag1255.PDS'; 
                    plxfile = '20191218/leo20191218_l3p2_final/leo20191218_l3p2.pl2';
                    ppNum = 35;
                    
            end
        otherwise
            error('falied to pick valid subect');
    end
    
    
    %% LOAD EXPERIMENT:
    if isempty(pdsfile)
        warning('NO SUCH PDS FILE EXISTS');
    end
    if isempty(plxfile)
        warning('NO SUCH PLX FILE EXISTS');
    end
    
    % load pds struct
    load([dirDataPDS pdsfile], '-mat');
    plxfile = [dirDataPlx plxfile];
    % get necessary ephys timing and timing conversions
    %    [PL2PTBfit,PL2PTB,PTB2PL, maxreconstructionerror, flagPlexonTrialStartTime] = ephys.syncPlexonClock(PDS, plxfile);
    syncClock = ephys.syncPlexonClock(PDS, plxfile);
    
end