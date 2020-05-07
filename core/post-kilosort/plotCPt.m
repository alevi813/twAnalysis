function [cpm, cpse, psthTime] = plotCPt(S, dataPath, spNum)
% calculate and plot avg CPt from a given group of cells (S)

if nargin <3
    spNum = [];
end

%%

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
        dirprob        = stim.dirprob(stim.validTrials);
        targchosen     = stim.targchosen(stim.validTrials);
        trialId        = stim.trialId(stim.validTrials);
    else
        dirprob        = stim.dirprob(stim.goodtrial);
        targchosen     = stim.targchosen(stim.goodtrial);
        trialId        = stim.trialId(stim.goodtrial);
    end
    
    % get a revco index, but only use it for traditional_cp unless useAll
    % flag is false.
    isRevco    = dirprob==0;
    isFrozen   = trialId==stim.frozenTrialIds;
    
    spikes_true = S(iNeuron).model(1).trialSpikes;
    rates_true  = S(iNeuron).model(1).trialRates;
    
    psthTime = S(iNeuron).model(2).psthTime;
    
    %% check from nans in the trial spikes
    nanix = isnan(mean(spikes_true, 2));
    
    spikes_true(nanix, :) = [];
    rates_true(nanix,:)   = [];
    
    % save every for unit...
    alldat.spkTrue{iNeuron}  = spikes_true;
    alldat.rateTrue{iNeuron} = rates_true;
    
    
    %% calculate cp
    
    % find which target corresponded to positive pulses
    % should be 2
    targPlus = mode(stim.targcorrect(sum(sum(stim.pulses, 3),2)>0));
    
    targPref = targPlus;
    % if you have d' < 0, assign targPref = 1
    if S(iNeuron).model(1).dprime < 0
        targPref = setdiff([1 2], targPlus);
        
        targchosen(targchosen ==2) = 0;
    else
        targchosen(targchosen ==1) = 0;
        targchosen(targchosen ==2) = 1;
    end
    
    % clean up some NaNs
    targchosen(nanix, :)       = [];
    isRevco(nanix)             = [];
    dirprob(nanix)             = [];
    isFrozen(nanix)            = [];
    
    
    % reassign some things to be saved
    alldat.targchosen{iNeuron}       = targchosen;
    alldat.isRevco{iNeuron}          = isRevco;
    alldat.isFrozen{iNeuron}         = isFrozen;
    alldat.dirprob{iNeuron}          = dirprob;
    
    % cp
    if size(rates_true, 1) ~= length(isRevco)
        trialDiff = length(isRevco) -  size(rates_true, 1);
        
        isRevco  = isRevco(trialDiff+1:end);
        isFrozen = isFrozen(trialDiff+1:end);
        
        trMismatch(iNeuron) = true;
    else
        trMismatch(iNeuron) = false;
    end
    
        % traditional, all revco
        [traditional_cp.m(iNeuron, :), traditional_cp.s(iNeuron, :)] = choiceProbabilityCalculate(rates_true(isRevco,:), targchosen(isRevco));
        % traditional, frozen
        [traditional_cp.frozen.m(iNeuron, :), traditional_cp.frozen.s(iNeuron, :)] = choiceProbabilityCalculate(rates_true(isFrozen,:), targchosen(isFrozen));
    
    
%         % traditional, all revco
%         [traditional_cp.m(iNeuron, :), traditional_cp.s(iNeuron, :)] = choiceProbabilityCalculate(spikes_true(isRevco,:), targchosen(isRevco));
%         % traditional, frozen
%         [traditional_cp.frozen.m(iNeuron, :), traditional_cp.frozen.s(iNeuron, :)] = choiceProbabilityCalculate(spikes_true(isFrozen,:), targchosen(isFrozen));
    
%     pdur      = 150;
%     pulseTime = [-100 50 200 350 500 650 800 950 1100];
%     
%     for iPulse = 1:9
%         ptimeIx = psthTime>=pulseTime(iPulse) & psthTime<pulseTime(iPulse)+pdur;
%         pulseSpikes = spikes_true(:, ptimeIx);
%         tmpAUC = roc(sum(pulseSpikes(isFrozen, :)'), targchosen(isFrozen));
%         cpPulse.frozen.m(iNeuron, iPulse) = tmpAUC.AUC;
%         cpPulse.frozen.se(iNeuron, iPulse) = tmpAUC.serror;
%         %        [cpPulse.frozen.m(iNeuron, iPulse), cpPulse.frozen.s(iNeuron, iPulse)] = roc(pulseSpikes, targchosen(isFrozen));
%     end
end % cell loop

%%
% get some min/max info for the fancy rectangle
if size(traditional_cp.frozen.m, 1) > 1
    %cpm  = smooth(mean(traditional_cp.frozen.m), 15);
    cpm  = gaussianSmooth(mean(traditional_cp.frozen.m), 15);
    cpse = std(traditional_cp.frozen.m) / sqrt(size(traditional_cp.frozen.m, 1));
    cpse = gaussianSmooth(cpse, 15);
else
    cpm  = gaussianSmooth(traditional_cp.frozen.m, 15);
    cpse = gaussianSmooth(traditional_cp.frozen.s, 15);
end

mx = max(cpm+cpse);
mn = min(cpm-cpse);

toplim    = mx + 0.005;
botlim    = mn - 0.005;

% plot
if isempty(spNum)
    %figure
else
    subplot(3,3,spNum); hold on
end
% rectangle('Position', [0, mn-0.005, 1050, toplim- botlim], 'EdgeColor', [0.9 0.9 0.9], 'FaceColor', [0.9 0.9 0.9]);
plot([psthTime(1) psthTime(end)], [0.5 0.5], 'k--')


if strcmp(stim.condition, 'flat')
    %blues    
    clr = [0 0 0.6];    
elseif strcmp(stim.condition, 'late')
    %yellers
    clr = [0.6 0.6 0];
else strcmp(stim.condition, 'early')
    %reds
    clr = [.6 0 0];
end

%clr = [.4 .4 .4];

boundedline(psthTime, cpm(1:length(psthTime)), cpse(1:length(psthTime)), 'cmap', clr, 'alpha');

xlim([-500 1500])
%ylim([.45 .55])
axis square
set(gca, 'FontSize', 12)
set(gca, 'box', 'off')


