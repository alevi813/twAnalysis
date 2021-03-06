
% path stuff and neuron file names
condition = {'nancy'};
dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/' condition{1}];
nDataPath = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/' condition{1} '/neurons'];

nFiles = dir(nDataPath); %nFiles = nFiles(4:end);
nFiles=nFiles(~ismember({nFiles.name},{'.','..','.DS_Store'})); % clean up

nFilesNames = arrayfun(@(x) x.name(1:9) , nFiles, 'UniformOutput', false);

%%
% WHICH EXPERIMENT WOULD YOU LIKE TO LOAD?
exname = 'n20161013';

nNeurons = sum(strcmp(nFilesNames, exname));

for iN = 1:nNeurons
    if iN < 10
        neurons(iN) = load([nDataPath filesep exname '_0' num2str(iN)]);
    else
        neurons(iN) = load([nDataPath filesep exname '_' num2str(iN)]);
    end
end

% load stim file
stim    = getStim(exname, dataPath);

% stimulus timing variables, bin size, windows
motionOnset = [stim.timing(:).motionon] + [stim.timing(:).plxstart];
%goTime      = [stim.timing(:).fpoff] + [stim.timing(:).plxstart];
goTime      = [stim.timing(:).fpoff] + [stim.timing(:).plxstart]+ stim.sacTime';

binSize = 0.01; % 10ms (in seconds)
window  = [-.5 3]; % 500 ms before motion onset to ~1s after (in seconds)
%window  = [-2.5 0]; % backwards windown when aligned to 'go'

nBins = diff(window)/binSize;

figure
% loop over neurons in the session
for kNeuron = 1:nNeurons
    
    spikeTimes  = neurons(kNeuron).spikeTimes;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% get spikes aligned to MOTION ONSET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [motionOn.spcnt, motionOn.bins]  = pdsa.binSpTimes(spikeTimes, motionOnset, window, binSize);
    
    % smooth spike count with a boxcar filter
    sm = 5; % size of boxcar (number of bins)
    tmp = filter(boxcar(sm)/sm, 1, motionOn.spcnt');
    tmp = flipud(tmp);
    tmp = filter(boxcar(sm)/sm, 1, tmp);
    tmp = flipud(tmp);
    
    motionOn.spcnt = tmp';
    motionOn.sprate = motionOn.spcnt/binSize;
    %motionOn.sprate = motionOn.sprate(goodTrials, :);
    
    goodTrials  = stim.goodtrial & ~any(isnan(motionOn.spcnt), 2);
    froIx       = stim.trialId==stim.frozenTrialIds;
    froIx = froIx(goodTrials);
    
    motionOn.sprate = motionOn.sprate(goodTrials, :);
    %motionOn.sprate = motionOn.sprate(froIx, :);

    Cho = sign(stim.targchosen - 1.5);
    Cho = Cho(goodTrials);
    %Cho = Cho(froIx);

    %plot
    odds = 1:nNeurons*2; odds = odds(rem(odds, 2)==1);
%    subplot(nNeurons, 2, odds(kNeuron))
    subplot(2, nNeurons, kNeuron)
    plot(motionOn.bins, nanmean(motionOn.sprate(Cho==1, :) )); hold on
    plot(motionOn.bins, nanmean(motionOn.sprate(Cho==-1, :) ));
    %xlim([-.5 1.6])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% get spikes aligned to GO SIGNAL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [fromGo.spcnt, fromGo.bins]  = pdsa.binSpTimes(spikeTimes, goTime, [-1 0], binSize);
    
    % smooth spike count with a boxcar filter
    sm = 5; % size of boxcar (number of bins)
    tmp_fromGo = filter(boxcar(sm)/sm, 1, fromGo.spcnt');
    tmp_fromGo = flipud(tmp_fromGo);
    tmp_fromGo = filter(boxcar(sm)/sm, 1, tmp_fromGo);
    tmp_fromGo = flipud(tmp_fromGo);
    
    fromGo.spcnt = tmp_fromGo';    
    fromGo.sprate = fromGo.spcnt/binSize;
    fromGo.sprate = fromGo.sprate(goodTrials, :);
    
    %plot
    evens = 1:nNeurons*2; evens = evens(rem(evens, 2)==0);
%    subplot(nNeurons, 2, evens(kNeuron))
    subplot(2, nNeurons, kNeuron+nNeurons)
    plot(fromGo.bins, nanmean(fromGo.sprate(Cho==1,:))); hold on
    plot(fromGo.bins, nanmean(fromGo.sprate(Cho==-1,:)))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % neuron loop

