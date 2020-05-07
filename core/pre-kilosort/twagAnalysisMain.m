cd /Users/Aaron/Dropbox/twagAnalysis4.1
%% load each session in a cell 'ppAll'. Save each 'pp' file.
subject = 'nancy';
eList   = [81:103];

dirDataPP  = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/pp/';

[ppAll] = dataFactory(subject, eList);

for ii = 1:length(ppAll)
    [ppAll{ii}] = twagCompute(ppAll{ii});
    pp = ppAll{ii};
    supertitle(['session ' num2str(eList(ii))], 14)
    %% Save experiment
    fNamepp = [subject '-' sprintf('e%02.0f', eList(ii)) '-pp'];
    %save(fullfile(dirDataPP, fNamepp), 'pp');
    clear pp
end



%% combine datasets from multiple sessions
subject = 'nancy';
%eList   = [58 60 61 62 64 65 66 67 68 69]; %all flat w/ephys
%eList   = [70 71 72 74 75 76 77 78 79 80]; %all late w/ephys
%eList   = [70 71 72 74]; %late transition w/ephys
%eList   = [75 76 77 78 79 80]; %last 6 late w/ephys
%eList   = [64 65 66 67 68 69]; %last six FLAT w/ephys

eList   = [81:103];

[ppAll] = ppLoader(subject, eList);

[gg] = combineDatasets(ppAll);

clear cohPerPulse
clear pmfUnfolded
clear kernel

[gg] = twagCompute(gg);

% save gg
dirDataGG  = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/';
fNamepp = [subject '-gg'];
%save(fullfile(dirDataGG, fNamepp), 'gg');

%supertitle(['sessions ' num2str(eList(1)) ' through ' num2str(eList(length(eList)))], 14)