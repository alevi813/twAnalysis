clear all
cd /Users/Aaron/Dropbox/twagAnalysis4.1
%%
ggbl = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-BL.mat');

[k.b, ~, k.s] = computeKernel(ggbl.gg.cohmat, ggbl.gg.stimDistNum, ggbl.gg.targRchosen, 0);

figure
% fit and plot actual pmf from BASELINE data
pmf   = pmfTools(ggbl.gg.sumCoh, ggbl.gg.targRchosen, 'nBins', 7);
pmf   = fit(pmf);
%plot(pmf, 'Color', [1 0.1 0.1], 'LineWidth', 2,'MarkerSize', 25);
plot(pmf, 'Color', [0.7 0.5 0], 'LineWidth', 2, 'MarkerSize', 25);

cohmat = ggbl.gg.cohmat;
cohmat   = cohmat/std(cohmat(:));

% simulate bernoulli observer with BASELINE weights
 x0=[ones(size(cohmat,1),1) cohmat]*k.b;
%x0=[ones(size(cohmat,1),1) cohmat]* ones(8,1);

pright=1./(1+exp(x0));
wCho=rand(size(cohmat,1),1)>pright;

pmf = pmfTools(ggbl.gg.sumCoh, wCho', 'nBins', 7);
pmf = fit(pmf);
plot(pmf, 'Color', [0.6 0.6 0.6], 'LineWidth', 2,'MarkerSize', 25);

% format figure
xlim([-60 60])
set(gcf, 'color', 'w');
set(gcf, 'units', 'inches', 'pos', [0 0 4.5 5])

% simulate bernoulli observer with LATE weights
gglw = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-LW-6t14.mat');

[k.b, ~, k.s] = computeKernel(gglw.gg.cohmat, gglw.gg.stimDistNum, gglw.gg.targRchosen, 0);
x0=[ones(size(cohmat,1),1) cohmat]*k.b;

pright=1./(1+exp(x0));
wCho=rand(size(cohmat,1),1)>pright;

pmf = pmfTools(ggbl.gg.sumCoh, wCho', 'nBins', 7);
pmf = fit(pmf);
plot(pmf, 'Color', [0 0 0], 'LineWidth', 2, 'MarkerSize', 25);


% format figure
xlim([-60 60])
set(gcf, 'color', 'w');
set(gcf, 'units', 'inches', 'pos', [0 0 4.5 5])

%legend('Observed behavior', 'Simulated behavior - baseline weights', 'Simulated behavior - Late weights');

kernelProofFig = gcf;
set(kernelProofFig, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4.5 5]); 
figdir = '/Users/Aaron/Dropbox/twagAnalysis4.1/figures/';
saveas(kernelProofFig, fullfile(figdir, 'kernelProofBaselineStimulus3.pdf'))

%% do this again for the late biased cohmat
clear all
gglw = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-LW-6t14.mat');

[k.b, ~, k.s] = computeKernel(gglw.gg.cohmat, gglw.gg.stimDistNum, gglw.gg.targRchosen, 0);

figure
% fit and plot actual pmf from LATE data
pmf   = pmfTools(gglw.gg.sumCoh, gglw.gg.targRchosen, 'nBins', 7);
pmf   = fit(pmf);
%plot(pmf, 'Color', [1 0.1 0.1], 'LineWidth', 2, 'MarkerSize', 25);
plot(pmf, 'Color', [0.7 0.5 0], 'LineWidth', 2, 'MarkerSize', 25);

cohmat = gglw.gg.cohmat;
cohmat   = cohmat/std(cohmat(:));

% simulate bernoulli observer with LATE weights
x0=[ones(size(cohmat,1),1) cohmat]*k.b;

pright=1./(1+exp(x0));
wCho=rand(size(cohmat,1),1)>pright;

pmf = pmfTools(gglw.gg.sumCoh, wCho', 'nBins', 7);
pmf = fit(pmf);
plot(pmf, 'Color', [0 0 0], 'LineWidth', 2, 'MarkerSize', 25);

% simulate bernoulli observer with BASELINE weights
ggbl = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-BL.mat');

[k.b, ~, k.s] = computeKernel(ggbl.gg.cohmat, ggbl.gg.stimDistNum, ggbl.gg.targRchosen, 0);
x0=[ones(size(cohmat,1),1) cohmat]*k.b;

pright=1./(1+exp(x0));
wCho=rand(size(cohmat,1),1)>pright;

pmf = pmfTools(gglw.gg.sumCoh, wCho', 'nBins', 7);
pmf = fit(pmf);
plot(pmf, 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'MarkerSize', 25);


% format figure
set(gcf, 'color', 'w');
set(gcf, 'units', 'inches', 'pos', [0 0 4.5 5])
%legend('Observed behavior', 'Simulated behavior - baseline weights', 'Simulated behavior - Late weights');

kernelProofFig = gcf;
set(kernelProofFig, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 4.5 5]); 
figdir = '/Users/Aaron/Dropbox/twagAnalysis4.1/figures/';
saveas(kernelProofFig, fullfile(figdir, 'kernelProofLateStimulus3.pdf'))

%% compare the two models
% wp=k.b(2:end); wp=wp/sum(wp);
% 
% figure(10); clf
% pmf   = pmfTools(pp.cohmat*wp, pp.targRchosen(:), 'nBins', 7);
% pmf   = fit(pmf);
% plot(pmf, 'Color', 'r');
% 
% pmf   = pmfTools(mean(pp.cohmat,2), pp.targRchosen(:), 'nBins', 7);
% pmf   = fit(pmf);
% plot(pmf, 'Color', 'k');




