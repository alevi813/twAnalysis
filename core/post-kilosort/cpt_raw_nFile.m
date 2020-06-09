function [allCP, bins] = cpt_raw_nFile(condition)

experiments = getExperimentsAnd(condition);

comp = getComp;

if strcmp(comp, 'desktop')
    dataPath  = ['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/' condition{1}];
    nDataPath = ['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/' condition{1} '/neurons'];
    
else
    dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/' condition{1}];
    nDataPath = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/' condition{1} '/neurons'];
    
end

nFiles = dir(nDataPath); %nFiles = nFiles(4:end);
nFiles=nFiles(~ismember({nFiles.name},{'.','..','.DS_Store'})); % clean up

% get a list of all neuron file names
nFilesNames = arrayfun(@(x) x.name(1:9) , nFiles, 'UniformOutput', false);

%%
allCP.m = [];

for kEx = 1:numel(experiments)
    
    exname=experiments{kEx};
    %clear neurons
    % get neurons for the given session
    nNeurons = sum(strcmp(nFilesNames, exname));
    
    for iN = 1:nNeurons
        if iN < 10
            neurons(iN) = load([nDataPath filesep exname '_0' num2str(iN)]);
        else
            neurons(iN) = load([nDataPath filesep exname '_' num2str(iN)]);
        end
    end
    
    % load stim file
    stim      = getStim(exname, dataPath);
    goodtrial = stim.goodtrial;
    
    Cho   = sign(stim.targchosen - 1.5);
    Direc = sign(sum(sum(stim.pulses, 3),2));
    froIx = stim.trialId==stim.frozenTrialIds;

    if corr(Cho, Direc) < 0
        Cho = -Cho;
    end
    binaryCho = Cho==1;
    
    motionOnset  = [stim.timing(:).motionon] + [stim.timing(:).plxstart];
    motionOffset = [stim.timing(:).motionoff] + [stim.timing(:).plxstart];
    goTime       = [stim.timing(:).fpoff] + [stim.timing(:).plxstart];
    
    froIx  = froIx(goodtrial);
    goTime = goTime(goodtrial);
    motionOnset = motionOnset(goodtrial);
    binaryCho = binaryCho(goodtrial);
    
    binSize = 0.01; % 10ms (in seconds)
    window  = [-1 0];

     % loop over neurons in the session
    for kNeuron = 1:nNeurons
                
        spikeTimes   = neurons(kNeuron).spikeTimes;
        
        % get spike count aligned to motion onset
        %[spcnt, bins]  = pdsa.binSpTimes(spikeTimes, motionOnset, window, binSize);
        [spcnt, bins]  = pdsa.binSpTimes(spikeTimes, goTime, window, binSize);
        
        % smooth spike count with a boxcar filter
        sm = 5; % size of boxcar (number of bins)
        tmp = filter(boxcar(sm)/sm, 1, spcnt');
        tmp = flipud(tmp);
        tmp = filter(boxcar(sm)/sm, 1, tmp);
        tmp = flipud(tmp);
        
        spcnt = tmp';
        
        % conver to spike rate
        sprate = spcnt/binSize;

%         nanix = isnan(sprate(:,1));
%         
%         sprate(nanix,:) = [];
%         froIx(nanix)   = [];
%         binaryCho(binaryCho(nanix)) = [];
        
        [tmpCP] = choiceProbabilityCalculate(sprate(froIx, :), binaryCho(froIx));
        
        cp.m(kNeuron, :) = tmpCP;
    end % kNeuron
    
    allCP.m = [allCP.m; cp.m];
    clear cp
end % kEx

%boundedline(bins, nanmean(allCP.m), nanstd(allCP.m)/sqrt(size(allCP.m,1)), 'cmap', [0 0 .5], 'alpha');
