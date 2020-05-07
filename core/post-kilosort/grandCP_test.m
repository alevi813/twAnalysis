function cp = grandCP_test(S, dataPath, shuff)


if nargin < 3
   shuff = false; 
end

for iNeuron = 1:length(S)
    
    % get your trial spikes
%    trSpikes = S(iNeuron).model(1).trialSpikes;
    trSpikes = S(iNeuron).model(1).trialRates;
    
    % load stim file
    if strcmp(S(iNeuron).exname(1), 'n')
        stim = load([dataPath{1} filesep 'stim' filesep S(iNeuron).exname '_stim.mat']);
    else
        stim = load([dataPath{2} filesep 'stim' filesep S(iNeuron).exname '_stim.mat']);
    end
    
    % index of frozen trials
    froIx = stim.trialId==stim.frozenTrialIds;
    
    if isfield(stim, 'validTrials')
        froIx = froIx(stim.validTrials);
        cho   = stim.targchosen(stim.validTrials);
    else
        froIx = froIx(stim.goodtrial);
        cho   = stim.targchosen(stim.goodtrial);
    end
    
    % assign pref targ based on d'
    if S(iNeuron).model(1).dprime > 0
        cho   = cho - 1;
    else
        cho = abs(cho - 2);
    end
    
    % separate frozen trials before ROC
    froSpikes = trSpikes(froIx, 100:220);
    froCho    = cho(froIx);
    
    nFroTrials  = size(froSpikes, 1);
    nBins       = size(froSpikes, 2);
    
    if shuff
        shuffIx   = randperm(nFroTrials);
        tmpSpikes = nan(nFroTrials, nBins);
        for iT = 1:nFroTrials
            tmpSpikes(iT, :) = froSpikes(shuffIx(iT), :);
        end
        
        froSpikes = tmpSpikes;
    end
    
    % do ROC
    if isequal(size(trSpikes, 1), length(froIx))
        [~, ~, ~, auc_tmp] = perfcurve(froCho, mean(froSpikes, 2), 1);
        
        cp(iNeuron) = auc_tmp;
    else
        cp(iNeuron) = NaN;
    end
end