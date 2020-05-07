 function [combo] = aggBehaviorFromStim(condition)
% aggregate pmf and ppk from STIM files

%condition   = {'nancy', 'early'};
%condition   = {'early'};

experiments = getExperimentsAnd(condition);

comp = getComp;

if strcmp(comp, 'desktop')
    dataPath{1}  = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo/stim';
    dataPath{2}  = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/stim';
else
    dataPath{1}  = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo/stim';
    dataPath{2}  = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/nancy/stim';
end

% combin the vars you want across sessions
combo.cho       = [];
combo.pulses    = [];
combo.pulses_z  = [];
combo.froIx     = [];
combo.dirprob   = [];
combo.sumcoh    = [];
combo.sumcoh_z  = [];

combo.motionOn  = [];
combo.goTime    = [];

for iExp = 1:length(experiments)
    
    if strcmp(experiments{iExp}(1), 'l')
        stim = load([dataPath{1} filesep experiments{iExp} '_stim.mat']);
    else
        stim = load([dataPath{2} filesep experiments{iExp} '_stim.mat']);
    end
    
    pulses = sum(stim.pulses, 3) ./19;
    sumcoh = mean(pulses,2);
    froIx = stim.trialId==stim.frozenTrialIds;
    
    if isfield(stim, 'validTrials')
        combo.dirprob = [combo.dirprob; stim.dirprob(stim.validTrials)];
        combo.cho     = [combo.cho; stim.targchosen(stim.validTrials)];
        pulses        = pulses(stim.validTrials, :);
        sumcoh        = sumcoh(stim.validTrials, :);
        froIx         = froIx(stim.validTrials, :);
        cho           = stim.targchosen(stim.validTrials);
    else
        combo.dirprob = [combo.dirprob; stim.dirprob(stim.goodtrial)];
        combo.cho     = [combo.cho; stim.targchosen(stim.goodtrial)];
        pulses        = pulses(stim.goodtrial, :);
        sumcoh        = sumcoh(stim.goodtrial, :);
        froIx         = froIx(stim.goodtrial, :);
        cho           = stim.targchosen(stim.goodtrial);
    end
    
    motionon = arrayfun(@(x) x.motionon , stim.timing, 'UniformOutput', false);
    motionon = cell2mat(motionon);
    motionon = motionon(stim.goodtrial);

    goTime = arrayfun(@(x) x.fpoff , stim.timing, 'UniformOutput', false);  
    goTime = cell2mat(goTime);
    goTime = goTime(stim.goodtrial);

%     if strcmp(stim.exname(1), 'l')
%         pulses_z = pulses ./ 0.2946; % leo
%         sumcoh_z = sumcoh ./ 0.2946; % leo
%     else
%         pulses_z = pulses ./ 0.3669; % nancy        
%         sumcoh_z = sumcoh ./ 0.3669; % nancy
%     end
    
    combo.pulses    = [combo.pulses; pulses];
    combo.sumcoh    = [combo.sumcoh; sumcoh];
%     combo.pulses_z  = [combo.pulses_z; pulses_z];
%     combo.sumcoh_z  = [combo.sumcoh_z; sumcoh_z];
    combo.froIx     = [combo.froIx; froIx];
    
    combo.motionOn  = [combo.motionOn; motionon];
    combo.goTime    = [combo.goTime; goTime];
    
    %    figure
    ppk = ppkTools(zscore(pulses), cho-1, 'ridge', true);
    ppkSlope(iExp) = table2array(ppk.w_lm.Coefficients(2,1));
    %     plot(ppk, 'plotFit', false, 'color', [0.7 0.7 0]);
    %     axis square; clear cho; clear pulses;
end % iExp loop

combo.cho      = logical(combo.cho-1); % important! stim denotes targs as 1 and 2, but they need to be 0 and 1
combo.pulses_z = zscore(combo.pulses);
combo.sumcoh_z = zscore(combo.sumcoh);

%%

nbins = 12;
doNot = true;
if ~doNot
    % playing around with folded pmf
    pmf   = pmfTools(combo.sumcoh, combo.cho, 'nBins', nbins);
    
    y(5) = mean([1-pmf.y(1) pmf.y(10)]);
    y(4) = mean([1-pmf.y(2) pmf.y(9)]);
    y(3) = mean([1-pmf.y(3) pmf.y(8)]);
    y(2) = mean([1-pmf.y(4) pmf.y(7)]);
    y(1) = mean([1-pmf.y(5) pmf.y(6)]);
    
    tmpse = mean(pmf.err);
    se(5) = mean([tmpse(1) tmpse(10)]);
    se(4) = mean([tmpse(2) tmpse(9)]);
    se(3) = mean([tmpse(3) tmpse(8)]);
    se(2) = mean([tmpse(4) tmpse(7)]);
    se(1) = mean([tmpse(5) tmpse(6)]);
    
    errorbar(pmf.x(6:end), y, se, 'v', 'color', [.4 .4 .4]); hold on
    [xfit, fitvals] = logFit(pmf.x(6:end), y);
    plot(xfit, fitvals, '--', 'color', [.4 .4 .4], 'LineWidth', 1.25)
    axis square
end

%%
%figure(4);
%subplot(3,3,9)
subplot(2,3,4)
pmf   = pmfTools(combo.sumcoh_z, combo.cho, 'nBins', nbins);
pmf   = fit(pmf);
pmf   = plot(pmf, 'color', [.6 0 0]);
axis square
box off

%figure(2);
%subplot(3,3,9)
subplot(2,3,1)
ppk = ppkTools(combo.pulses_z, combo.cho);
plot(ppk, 'plotFit', false, 'color', [.6 0 0]);

hold on
ppk_rc = ppkTools(combo.pulses_z(combo.dirprob==0, :), combo.cho(combo.dirprob==0));
plot(ppk_rc, 'plotFit', false, 'color', [.6 0 0], 'linestyle', '--');

ylim([0 .8])
ylim([0 .7])
xticks([0 1 2 3 4 5 6 7 8])
axis square

combo.pmf = pmf;
combo.ppk = ppk;
%% slopes & 95%ci
% early: -0.10865  +/- 0.0275
% flat : -0.090801 +/- 0.0223
% late : 0.083163 +/- 0.0678