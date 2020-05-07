ggbl = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-BL-57t69.mat');

%%
winSize = 500;
nTrials = length(ggbl.gg.cohmat);
figure; hold on

cohmat      = ggbl.gg.cohmat;
stimDistNum = ggbl.gg.stimDistNum;
targRchosen = ggbl.gg.targRchosen;

% %%%% shuffle trials if you want %%%%
% nTr = length(cohmat);
% randInd = randperm(nTr);
% 
% cohmat      = cohmat(randInd, :);
% stimDistNum = stimDistNum(randInd);
% targRchosen = targRchosen(randInd);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:nTrials-winSize
    tempPulses      = cohmat(ii:(ii+winSize-1),:);
    tempStimDist    = stimDistNum(ii:(ii+winSize-1));
    tempTargRchosen = targRchosen(ii:(ii+winSize-1));
    
    [k, ~, ~] = computeKernel(tempPulses, tempStimDist, tempTargRchosen, 1);
    
    w = (sum(k(1:3)) - sum(k(5:7))) / (sum(k(1:3)) + sum(k(5:7)));
    
    %plot(ii, w, 'ko')
    plot(ii, w, 'o', 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6])
end

%% late biased
%clear all
jj = ii;

gglw = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-LATE-70t77.mat');
%gglw = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-LW-6t12.mat');
%gglw = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-LW-6t14.mat');
winSize = 500;
nTrials = length(gglw.gg.cohmat);


cohmat      = gglw.gg.cohmat;
stimDistNum = gglw.gg.stimDistNum;
targRchosen = gglw.gg.targRchosen;

% %%%% shuffle trials if you want %%%%
% nTr = length(cohmat);
% randInd = randperm(nTr);
% 
% cohmat      = cohmat(randInd, :);
% stimDistNum = stimDistNum(randInd);
% targRchosen = targRchosen(randInd);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:nTrials-winSize
    tempPulses      = cohmat(ii:(ii+winSize-1),:);
    tempStimDist    = stimDistNum(ii:(ii+winSize-1));
    tempTargRchosen = targRchosen(ii:(ii+winSize-1));
    
    [k, ~, ~] = computeKernel(tempPulses, tempStimDist, tempTargRchosen, 1);
    
    w = (sum(k(1:3)) - sum(k(5:7))) / (sum(k(1:3)) + sum(k(5:7)));
    
    %plot(jj, w, 'ro')
    plot(jj, w, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0])

    jj = jj+1;
end

ylabel('Weighting Index', 'FontSize', 16)
set(gca, 'FontSize', 14)
xlabel('Trial window start', 'FontSize', 16) 
set(gca, 'FontSize', 14)

hh = refline(0,0);
set(hh, 'color', 'k', 'linestyle', '--');

wSessionsFig = gcf;
set(wSessionsFig, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 3.7]); 
figdir = '/Users/Aaron/Dropbox/twagAnalysis4.1/';
saveas(wSessionsFig, fullfile(figdir, 'wSessionsFig.pdf'))






%% trying out a way to do combined jake + aaron data for nancy

% get and combine data 
clear all

cd('/Users/Aaron/Dropbox/PatStimForAaron');
setupPaths

%import baseline data
figDir='./figures/behavior';
import behavior.*
import regression.*
S = BehaviorSummary();

% baseline

X=cell2mat(arrayfun(@(x) x.X', S, 'Uniformoutput', false))';
cho=cell2mat(arrayfun(@(x) x.Y', S, 'Uniformoutput', false))'==1;
coh=cell2mat(arrayfun(@(x) x.coh', S, 'Uniformoutput', false))';
trix=cell2mat(arrayfun(@(x) x.isrevco', S, 'Uniformoutput', false))'==1;
nix=cell2mat(arrayfun(@(x) x.isNancy.*true(1,numel(x.isrevco)), S, 'Uniformoutput', false))'==1;

% baseline data from Jake
ix1=nix(:);
ix0=trix(:)&nix(:);

% baseline data from Aaron
load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-BL.mat')
cohmat = gg.cohmat(gg.stimDistNum==3, :);
cohmat = cohmat/std(cohmat(:));
targCho = gg.targRchosen(gg.stimDistNum==3)';

% combine BL data
bigCohmat = [X(ix0,:); cohmat];
bigCho    = [cho(ix0); targCho];

%%
nTrials = length(bigCohmat);
winSize = 500;
figure; hold on

for ii = 1:nTrials-winSize
    tempPulses      = bigCohmat(ii:(ii+winSize-1),:);
    %tempStimDist    = stimDistNum(ii:(ii+winSize-1));
    tempTargRchosen = bigCho(ii:(ii+winSize-1));
    
    %[k, ~, ~] = computeKernel(tempPulses, tempStimDist, tempTargRchosen, 1);
    
    [k,dev,s]=glmfit(tempPulses,tempTargRchosen,'binomial');
    
    w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));
    
    %plot(ii, w, 'ko')
    plot(ii, w, 'o', 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6])
end