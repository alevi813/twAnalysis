function [xAll, allAuc, aucse] = nmfByNeuron(S, dataPath)

for iNeuron = 1:length(S)
    
    % load stim
    % get stim file from appropriate monkey path
    if strcmp(S(iNeuron).exname(1), 'n')
        stim = load([dataPath{1} filesep 'stim' filesep S(iNeuron).exname '_stim.mat']);
    else
        stim = load([dataPath{2} filesep 'stim' filesep S(iNeuron).exname '_stim.mat']);
    end
    
    condition{iNeuron} = stim.condition;
    
    if isfield(stim, 'validTrials')
        pulses = stim.pulses(stim.validTrials,:,:);
    else
        pulses = stim.pulses(stim.goodtrial,:,:);
    end
    
    pulses = sum(pulses, 3) ./19;
    sumcoh = mean(pulses,2);
    
    nBins  = 10;
    cohbins = binCoherences(sumcoh, nBins); % bin cohs for later
    
    direction = sign(sumcoh);
    noDir     = direction==0;
    direction = direction(~noDir) ==1;
    sumcoh    = sumcoh(~noDir);
    binId     = cohbins.id(~noDir);
    
    spikes = nanmean(S(iNeuron).model(1).trialRates(~noDir, (100:220)), 2);
    
    binsub = [0 1 2 3 4];
    
    for kBin = 1:nBins/2
        tmpdir    = [direction(binId==kBin); direction(binId==(nBins-binsub(kBin)))];
        tmpspikes = [spikes(binId==kBin); spikes(binId==(nBins-binsub(kBin)))];
        
        if length(unique(tmpdir)) > 1
            %[~, ~, ~, bin_auc(kBin)] = perfcurve(tmpdir, tmpspikes, 1);
            tmp(kBin) = roc(tmpspikes, tmpdir);
            
            bin_auc(kBin)   = tmp(kBin).AUC;
            
            if length(S) ==1
                aucse(kBin) = tmp(kBin).serror;
            else
                aucse(kBin) = NaN;                
            end
        else
            bin_auc(kBin) = NaN;
        end
        
        if bin_auc(kBin) < .5
            bin_auc(kBin) = 1-bin_auc(kBin);
        end
    end %kBin
    
    % fold/unsign your cohbins
    fitvals = abs(cohbins.binCenters(1:nBins/2));
    % save single expt vals into larger matrices
    xAll(iNeuron, :) = fitvals;
    allAuc(iNeuron, :) = bin_auc;
    
end %iNeuron