
%exDir = '/Users/Aaron/Desktop/test-nFile';
baseDir = '/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/nancy';

sessions = dir(baseDir);
sessions=sessions(~ismember({sessions.name},{'.','..','.DS_Store','allSessions.mat'})); % CLEAN DAT SHIT UP, also keep it from crashing if they already ran it and have a sessions.mat

for iSession = 9
%for iSession = 6
   
    baseDir = '/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/nancy';
    
    sessions = dir(baseDir);
    sessions=sessions(~ismember({sessions.name},{'.','..','.DS_Store','allSessions.mat'})); % CLEAN DAT SHIT UP, also keep it from crashing if they already ran it and have a sessions.mat
    
    %% preprocess a bunch of stuff
    dateDir = ([baseDir filesep sessions(iSession).name]);
    expts = dir(dateDir);
    expts = expts(~ismember({expts.name},{'.','..','.DS_Store'})); % CLEAN UP
    
    % get 'final' version of sort
    for ii = 1:length(expts)
        isFinal(ii) = strcmp(expts(ii).name(end-4:end), 'final');
    end
    
    expts = expts(isFinal);
    
    exDir = ([dateDir filesep expts(1).name]);
    
    load([exDir filesep 'PDS.mat']);
    load([exDir filesep 'syncClock.mat']);
    load([exDir filesep 'KiloSort' filesep 'rez.mat']);
    
    % set paths to spike times and clusters, load em.
    npy.spikeTimes    = ([exDir filesep 'KiloSort' filesep 'spike_times.npy']);
    npy.spikeClusters = ([exDir filesep 'KiloSort' filesep 'spike_clusters.npy']);
    
    % get spike times and clusters
    sampleRate = rez.ops.maxFR;
    st = readNPY(npy.spikeTimes);
    st = double(st);
    st = st./sampleRate;
    
    clusters   =  readNPY(npy.spikeClusters);
    
    % read in cluster info from the csv file with group labels, if it exists.
    if exist([exDir filesep 'KiloSort' filesep 'cluster_groups.csv'], 'file')
        [clusterID, clusterGroup] = readClusterCSV(exDir);
        clusterGood = clusterID(strcmp(clusterGroup, 'good'));
        clusterMua  = clusterID(strcmp(clusterGroup, 'mua'));
        
        clusterID = sort([clusterGood; clusterMua]);
    else
        clusterID  = unique(clusters);
    end
    nClust = length(clusterID);
    
    % get necessary vars from PDS
    goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
    goodtrial = cell2mat(goodtrial);
    goodtrial(isnan(goodtrial)) = 0;
    goodtrial = logical(goodtrial);
    trialnum = 1:length(goodtrial);
    
    % set the centering time based on motion start
    motionStart = cellfun(@(x) x.stimulus.statesStartTime(4),  PDS.data, 'UniformOutput', false);
    motionStart = cell2mat(motionStart)';
    
    if size(nsync.computerTrialStartTime) == size(motionStart)
        centerTime = nsync.plexonTrialStartTime + motionStart;
    else
        centerTime = nsync.plexonTrialStartTime + motionStart(1:length(nsync.computerTrialStartTime), :);
    end
    
    %% ok, now set up data in jake's format
    for iClust = 1:nClust
        n.brainArea = 'MT';
        %channel =
        n.exname = [PDS.initialParametersMerged.session.file(1) PDS.initialParametersMerged.session.file(6:13)];

        n.spikeTimes = double(st(clusters == clusterID(iClust)));
        for iTrial = 1:length(goodtrial)
            n.spikeCount(iTrial) = sum(n.spikeTimes > centerTime(iTrial) & n.spikeTimes < (centerTime(iTrial)+1.1));
        end
        
        % trial index logical
        n.trialIndex = goodtrial & n.spikeCount > 0;
        
        fullTrialIndex(:, iClust) = n.trialIndex;
    end % cluster loop
    figure
    imagesc(fullTrialIndex)
    clear all
end % session loop