function [pp] = twagCompute(pp, plotOn)
%% twagCompute
% description of function

%% cohPerPulse stats
[pp.cohPerPulse] = cohPerPulseStats(pp.cohmat, pp.stimDistNum);

%% compute and plot PMF

%pp.pmfUnfolded   = pmfTools(zscore(pp.sumCoh), pp.targRchosen, 'nBins', 18);
pp.pmfUnfolded   = pmfTools(pp.sumCoh, pp.targRchosen, 'nBins', 18);
pp.pmfUnfolded   = fit(pp.pmfUnfolded);

%% compute psychophysical kernel
% calculate the kernel. Set last arg = 1 to normalize by L2 norm.
% [pp.kernel.b, ~, pp.kernel.s] = computeKernel(pp.cohmat, pp.stimDistNum, pp.targRchosen, 1);

pp.ppk = ppkTools(pp.cohmat, pp.targRchosen', 'ridge', true);

%% plot stimulus statistics summary
if plotOn
figName = 'StimulusStatistics';
hFig    = figure('Name', figName, 'Position', [200 300 1000 1000], 'PaperSize', [12 16], 'PaperPosition', [0 0 12 16]);

% Coherence (full trial, aka sumCoh) distribution of overall session
subplot(221); hold on
plot.twagCohDist(pp.sumCoh, pp.stimDistNum, 'overall')

% Coherence (full trial, aka sumCoh) by 5 separate distributions
subplot(223); hold on
plot.twagCohDist(pp.sumCoh, pp.stimDistNum, 'byDist')

% Coherence per pulse for overall session
subplot(222); hold on
plot.visCohPerPulse(pp.cohPerPulse, 'overall')

% Coherence per pulse by distribution
subplot(224); hold on
plot.visCohPerPulse(pp.cohPerPulse, 'byDist')

% supertitle(pp.initialParametersMerged.session.file(1:11), 13);

%% plot behavior summary
figName = 'Behavior';
hFig    = figure('Name', figName, 'Position', [500 300 700 350], 'PaperSize', [12 16], 'PaperPosition', [0 0 10 5]);

% plot unfolded PMF:
subplot(121); hold on
pp.pmfUnfolded   = plot(pp.pmfUnfolded);
xlim([-60 60])
%title('Unfolded PMF', 'FontSize', 14);
xlabel('% Coherence', 'FontSize', 16);
set(gca, 'FontSize', 14)
ylabel('Right choices', 'FontSize', 16)
set(gca, 'FontSize', 14)

% plot kernel
subplot(122); hold on
% errorbar(1:7, pp.kernel.b, pp.kernel.s.se, ...
%     'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
% text(6, max(pp.kernel.b)+mean(pp.kernel.s.se), ['n = ' num2str(length(pp.cohmat(pp.stimDistNum == 3)))])
plot(pp.ppk, 'plotFit', false);

% make it pretty
xlabel('Pulse', 'FontSize', 16)
set(gca, 'FontSize', 14)
ylabel('Weight (normalized)', 'FontSize', 16)
set(gca, 'Xtick', 1:7)
set(gca, 'FontSize', 14)
xlim([.5 7.5]);

%supertitle(pp.initialParametersMerged.session.file(1:11), 13);
end
