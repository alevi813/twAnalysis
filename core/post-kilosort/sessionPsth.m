%% sessionPSTH --- to pair with mtrfmapPlx output (tuning curves)
comp = getComp;

% load in a thIx from mtrfmapPlx
%thIx = [13 12 11 14 16 4 6 17 5 7 8 9 2 1 10 3 15]; % s5_1
%thIx = [10 6 2 1 15 4 13 7 9 12 16 17 3 14 11 5 8];  % s5_2

%thIx = [3 11 1 5 12 16 4 8 6 18 14 9 13 15 17 2 7 10]; % s28_1
thIx = [9 4 8 16 6 10 18 14 15 17 2 13 1 7 3 5 12 11]; % s28_2

subject = 'leo';
if strcmp(subject, 'leo')
    fullDir = ('/Volumes/HukLab/Macaque/Projects/twag/leo');
else
    fullDir = ('/Volumes/HukLab/Macaque/Projects/twag/nancy');
end

sessions = dir(fullDir);
sessions=sessions(~ismember({sessions.name},{'.','..','.DS_Store','allSessions.mat'})); % CLEAN DAT SHIT UP, also keep it from crashing if they already ran it and have a sessions.mat

sessionList = 28; % be careful about sessions in your directory with no twagPDS file.
sessions = sessions(sessionList);

%screenLatency = .0524634; %52 ms screen latency for old LG/rig1 set up.
screenLatency = 0;


dateDir = ([fullDir filesep sessions.name]);
expts = dir(dateDir);
expts = expts(~ismember({expts.name},{'.','..','.DS_Store'})); % CLEAN UP

%     % there are one or two sessions with more than one expt. Currently hardcoding to only look at first right now.
%     if length(expts) > 1
%         expts = expts(1);
%     end

for ii = 1:length(expts)
    isFinal(ii) = strcmp(expts(ii).name(end-4:end), 'final');
    %isFinal(ii) = strcmp(expts(ii).name(end-3:end), 'test');
end

if size(expts, 1) > 1
    expts = expts(isFinal);
end

%% load  stuff
%  PDS, clock info, & rez
exDir = ([dateDir filesep expts.name]);

if strcmp(comp, 'desktop')
    stim = load(['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo/stim/' [subject(1) dateDir(end-7:end)] '_stim.mat']);
else
    stim = load(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo/stim/' [subject(1) dateDir(end-7:end)] '_stim.mat']);
end

% annoying workaround for loading pp...
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
load([exDir filesep f.name], '-mat'); % load PDS

% process PDS
pp = analyzeBehaviorSession(exDir, f.name, 1, false);


%% get spike times and clusters from n-files (formerly from kilosort directory)

if strcmp(comp, 'laptop')
    neuronDir       = dir('/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo/neurons');
else
    neuronDir       = dir('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo/neurons');
end
neuronFileNames = arrayfun(@(x) x.name, neuronDir, 'UniformOutput', false);
neuronDir       = neuronDir(contains(neuronFileNames, dateDir(end-7:end), 'IgnoreCase', true));

nClust = size(neuronDir, 1);

% get necessary vars from PDS
%         goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
%         goodtrial = cell2mat(goodtrial);
%         goodtrial(isnan(goodtrial)) = 0;
goodtrial = pp.goodtrial;

stimDistNum = cellfun(@(x) x.theseGabs.stimDistNum, PDS.data, 'UniformOutput', false);
stimDistNum = cell2mat(stimDistNum);

sumCoh = cellfun(@(x) x.theseGabs.sumCoh, PDS.data, 'UniformOutput', false);
sumCoh = cell2mat(sumCoh);

motionOn   = cell2mat(arrayfun(@(x) x.motionon, stim.timing, 'UniformOutput', false));
plxTrStart = cell2mat(arrayfun(@(x) x.plxstart, stim.timing, 'UniformOutput', false));

centerTime = plxTrStart+motionOn;

%% psth by motion distribution (let's make this its own function!!!)
figure2
nrows = 3;
ixix = [1 2 3 4 5];
clr  = [1 0 0; 1 .6 .6; 0 0 0; .6 .6 1; 0 0 1];
for iClust = 1:nClust
    %figure
    %theseSpikes = double(spikeTimesPtb(clusters == clusterID(thIx(iClust))));
    thisNeuron = load([neuronDir(iClust).folder filesep neuronDir(thIx(iClust)).name]);
    theseSpikes = thisNeuron.spikeTimes;

    maxrate(iClust) = 0 ; % preset max rate to be filled at end of dist plot loop
    
    for kk = 1:length(ixix)
        distIx   = stimDistNum(:)==ixix(kk) & goodtrial(:) == 1;
        if length(distIx) ~= length(centerTime)     % fix for sessions that run longer after turning off recording
            distIx = distIx(1:length(centerTime));
        end
        
        if sum(theseSpikes) > 0 % whhhhhyyyyyyyyY??
            [m,s,bc,v,~] = pdsa.eventPsth(theseSpikes, centerTime(distIx), [-.5 1.5],  .001, ones(100,1)/100);
        end
        
        hold on
        subplot(nrows, 6, iClust);
        plot(bc, m, 'Color', clr(kk,:), 'LineWidth', 1.75);
        xlim([-.25 1.5])
        %supertitle(num2str(dateDir(end-7:end)), 12);
        
        if maxrate(iClust) < max(m)
            maxrate(iClust) = max(m);
        end
    end
    %title(['ID ' num2str(clusterID(ii))]);
    box off        
end % cluster loop for plotting by dist
