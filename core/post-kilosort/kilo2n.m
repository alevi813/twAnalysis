
%exDir = '/Users/Aaron/Desktop/test-nFile';
%baseDir = '/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/leo';
baseDir = '/Volumes/HukLab/Macaque/Projects/twag/leo';

sessions = dir(baseDir);
sessions=sessions(~ismember({sessions.name},{'.','..','.DS_Store','allSessions.mat', '20191024', '20191121'})); % CLEAN DAT SHIT UP, also keep it from crashing if they already ran it and have a sessions.mat

%for iSession = 1:36
for iSession = 34:35
    
    baseDir = '/Volumes/HukLab/Macaque/Projects/twag/leo';
    
    sessions = dir(baseDir);
    sessions=sessions(~ismember({sessions.name},{'.','..','.DS_Store','allSessions.mat', 'test', 'toSort', '20191024', '20191121'})); % CLEAN DAT SHIT UP, also keep it from crashing if they already ran it and have a sessions.mat
    
    %% preprocess a bunch of stuff
    dateDir = ([baseDir filesep sessions(iSession).name]);
    expts = dir(dateDir);
    expts = expts(~ismember({expts.name},{'.','..','.DS_Store'})); % CLEAN UP
    
    % get 'final' version of sort
    if length(expts) > 1
        for ii = 1:length(expts)
            isFinal(ii) = strcmp(expts(ii).name(end-4:end), 'final');
        end
        
        expts = expts(isFinal);
    end
    
    exDir = ([dateDir filesep expts(1).name]);
    
    %load([exDir filesep 'PDS.mat']);
    load([exDir filesep 'syncClock.mat']);
    load([exDir filesep 'KiloSort' filesep 'rez.mat']);
    
    % annoying workaround for loading pp/PDS...
    f = dir(exDir);
    f = f(~ismember({f.name},{'.','..','.DS_Store'})); % CLEAN UP
    
    for ii = 1:length(f)
        f(ii).isPDS = strcmp(f(ii).name(end-3:end), '.PDS');
    end
    isPDS = arrayfun(@(x) x.isPDS, f, 'uniformoutput', 0);
    isPDS = cell2mat(isPDS);
    f = f(isPDS);
    
    % if you have multiple pds files, use the largest one
    for ii = 1:length(f)
        byteSize(ii) = f(ii).bytes;
    end
    [~, sizeIx] = max(byteSize);
    f = f(sizeIx);
    
    %load([f.folder filesep f.name], '-mat'); % load PDS
    load([exDir filesep f.name], '-mat');
    
    % set paths to spike times and clusters, load em.
    npy.spikeTimes    = ([exDir filesep 'KiloSort' filesep 'spike_times.npy']);
    npy.spikeClusters = ([exDir filesep 'KiloSort' filesep 'spike_clusters.npy']);
    
    % get spike times and clusters
    % sampleRate = rez.ops.maxFR;
    sampleRate = rez.ops.fs;
    st = readNPY(npy.spikeTimes);
    st = double(st);
    st = st./sampleRate;
    
    clusters   =  readNPY(npy.spikeClusters);
    
    % read in cluster info from the csv file with group labels, if it exists.
    %     if exist([exDir filesep 'KiloSort' filesep 'cluster_groups.csv'], 'file')
    %         [clusterID, clusterGroup] = readClusterCSV(exDir);
    %         clusterGood = clusterID(strcmp(clusterGroup, 'good'));
    %         clusterMua  = clusterID(strcmp(clusterGroup, 'mua'));
    %
    %         clusterID = sort([clusterGood; clusterMua]);
    %     else
    clusterID  = unique(clusters);
    %    end
    nClust = length(clusterID);
    
    % get channels
    spikeTemplates = rez.st3(:,2);
    amplitudes = rez.st3(:,3);
    
    temps = permute(rez.Wraw,[3 2 1]); % looks like it was there, but in a different order, so have to rearrange
    
    channel = nan(1,length(clusterID));
    for c = 1:length(clusterID)
        clust     = clusterID(c);
        I         = clusters == clust;
        templates =  unique(spikeTemplates(I));
        
        % temps = [nTemplates, nTimePoints, nTempChannels]
        t = squeeze(range((temps(templates,:,:)),1)); % what dimension to range over???????
        
        maxWFperCh = max(t);
        maxWF = max(max(t));
        if ~all(maxWFperCh==0)
            peakCh = find(maxWFperCh==maxWF);
        else
            peakCh = NaN;
        end
        
        channel(c) = peakCh;
    end
    
    % get necessary vars from PDS
    goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
    goodtrial = cell2mat(goodtrial);
    goodtrial(isnan(goodtrial)) = 0;
    goodtrial = logical(goodtrial);
    trialnum = 1:length(goodtrial);
    
    % set the centering time based on motion start
    motionStart = cellfun(@(x) x.stimulus.statesStartTime(4),  PDS.data, 'UniformOutput', false);
    motionStart = cell2mat(motionStart)';
    
    if size(syncClock.computerTrialStartTime) == size(motionStart)
        centerTime = syncClock.plexonTrialStartTime + motionStart;
    else
        centerTime = syncClock.plexonTrialStartTime + motionStart(1:length(syncClock.computerTrialStartTime), :);
    end
    
    %    dat = load(['/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData/leoSession' num2str(iSession)]);
    dat = load(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData/leoSession' num2str(iSession)]);
    
    %% ok, now set up data in jake's format
    for iClust = 1:nClust
        n.brainArea = 'MT';
        n.channel   = channel(iClust);
        n.id        = iClust;
        n.exname    = [PDS.initialParametersMerged.session.file(1) PDS.initialParametersMerged.session.file(4:11)];        
        
        n.spikeTimes = double(st(clusters == clusterID(iClust)));
        %        for iTrial = 1:length(goodtrial)
        for iTrial = 1:length(centerTime)
            n.spikeCount(iTrial) = sum(n.spikeTimes > centerTime(iTrial) & n.spikeTimes < (centerTime(iTrial)+1.1));
        end
        
        % trial index logical
        gtx          = goodtrial(1:length(centerTime));
        n.trialIndex = gtx & n.spikeCount > 0;
        
        % use trial index logical to get what you really want.
        n.spikeCount = n.spikeCount(n.trialIndex);
        n.trialIndex = trialnum(n.trialIndex);
        
        n.goodClustIx = dat.goodClustIx(iClust);
        %        save(['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo/neurons' filesep n.exname '_0' num2str(iClust)], '-struct', 'n', '-v7.3')
        save(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo/neurons' filesep n.exname '_0' num2str(iClust)], '-struct', 'n', '-v7.3')
    end % cluster loop
    clear all
end % session loop