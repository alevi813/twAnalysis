% get and combine data from jake and aaron
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
winSize = 500;
numBlocks = floor(length(bigCho)/winSize);

blockCoh = bigCohmat(1:winSize, :);
blockCho = bigCho(1:winSize);

[k,~,~]=glmfit(blockCoh,blockCho,'binomial');

w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));

figure; hold on
plot(1, w, 'o', 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6]);

for ii = 1:numBlocks-1
    blockCoh = bigCohmat((winSize*ii+1):(winSize*ii+winSize), :);
    blockCho = bigCho((winSize*ii+1):(winSize*ii+winSize));
    
    [k,~,~]=glmfit(blockCoh,blockCho,'binomial');
    
    w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));
    
    plot(ii+1, w, 'o', 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6]);
    
end

blockCoh = bigCohmat((winSize*ii+1):end, :);
blockCho = bigCho((winSize*ii+1):end);
[k,~,~]=glmfit(blockCoh,blockCho,'binomial');

w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));

plot(ii+2, w, 'o', 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6]);

%%
cd('/Users/Aaron/Dropbox/twagAnalysis4.1');
clear all

jj = 8;
gglw = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-6t12-L.mat');
%gglw = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-LW-6t14.mat');
cohmat      = gglw.gg.cohmat;
stimDistNum  = gglw.gg.stimDistNum;
targRchosen = gglw.gg.targRchosen;

winSize = 500;
numBlocks = floor(length(targRchosen)/winSize);

blockCoh     = cohmat(1:winSize, :);
blockCho     = targRchosen(1:winSize);
blockStimDistNum  = stimDistNum(1:winSize);

[k, ~, ~] = computeKernel(blockCoh, blockStimDistNum, blockCho, 1);

w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));

plot(jj, w, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);

for ii = 1:numBlocks-1
    blockCoh    = cohmat((winSize*ii+1):(winSize*ii+winSize), :);
    blockCho    = targRchosen((winSize*ii+1):(winSize*ii+winSize));
    blockStimDistNum = stimDistNum((winSize*ii+1):(winSize*ii+winSize));
    
    [k, ~, ~] = computeKernel(blockCoh, blockStimDistNum, blockCho, 1);
    
    w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));
    
    plot(jj+1, w, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
   
    jj= jj+1;
end

blockCoh     = cohmat((winSize*ii+1):end,:);
blockCho     = targRchosen((winSize*ii+1):end);
blockStimDistNum  = stimDistNum((winSize*ii+1):end);

[k, ~, ~] = computeKernel(blockCoh, blockStimDistNum, blockCho, 1);

w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));

plot(jj+1, w, 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);

%%
cd('/Users/Aaron/Dropbox/twagAnalysis4.1');
clear all

jj = 18;
gg = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-15t21.mat');

cohmat      = gg.gg.cohmat;
stimDistNum  = gg.gg.stimDistNum;
targRchosen = gg.gg.targRchosen;

winSize = 500;
numBlocks = floor(length(targRchosen)/winSize);

blockCoh     = cohmat(1:winSize, :);
blockCho     = targRchosen(1:winSize);
blockStimDistNum  = stimDistNum(1:winSize);

[k, ~, ~] = computeKernel(blockCoh, blockStimDistNum, blockCho, 1);

w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));

plot(jj, w, 'o', 'Color', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]);

for ii = 1:numBlocks-1
    blockCoh    = cohmat((winSize*ii+1):(winSize*ii+winSize), :);
    blockCho    = targRchosen((winSize*ii+1):(winSize*ii+winSize));
    blockStimDistNum = stimDistNum((winSize*ii+1):(winSize*ii+winSize));
    
    [k, ~, ~] = computeKernel(blockCoh, blockStimDistNum, blockCho, 1);
    
    w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));
    
    plot(jj+1, w, 'o', 'Color', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]);
   
    jj= jj+1;
end

blockCoh     = cohmat((winSize*ii+1):end,:);
blockCho     = targRchosen((winSize*ii+1):end);
blockStimDistNum  = stimDistNum((winSize*ii+1):end);

[k, ~, ~] = computeKernel(blockCoh, blockStimDistNum, blockCho, 1);

w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));

plot(jj+1, w, 'o', 'Color', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]);

%%
cd('/Users/Aaron/Dropbox/twagAnalysis4.1');
clear all

jj = 34;
gg = load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-22t57.mat');

cohmat      = gg.gg.cohmat;
stimDistNum  = gg.gg.stimDistNum;
targRchosen = gg.gg.targRchosen;

winSize = 500;
numBlocks = floor(length(targRchosen)/winSize);

blockCoh     = cohmat(1:winSize, :);
blockCho     = targRchosen(1:winSize);
blockStimDistNum  = stimDistNum(1:winSize);

[k, ~, ~] = computeKernel(blockCoh, blockStimDistNum, blockCho, 1);

w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));

plot(jj, w, 'o', 'Color', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]);

for ii = 1:numBlocks-1
    blockCoh    = cohmat((winSize*ii+1):(winSize*ii+winSize), :);
    blockCho    = targRchosen((winSize*ii+1):(winSize*ii+winSize));
    blockStimDistNum = stimDistNum((winSize*ii+1):(winSize*ii+winSize));
    
    [k, ~, ~] = computeKernel(blockCoh, blockStimDistNum, blockCho, 1);
    
    w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));
    
    plot(jj+1, w, 'o', 'Color', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]);
   
    jj= jj+1;
end

blockCoh     = cohmat((winSize*ii+1):end,:);
blockCho     = targRchosen((winSize*ii+1):end);
blockStimDistNum  = stimDistNum((winSize*ii+1):end);

[k, ~, ~] = computeKernel(blockCoh, blockStimDistNum, blockCho, 1);

w = (sum(k(2:4)) - sum(k(6:8))) / (sum(k(2:4)) + sum(k(6:8)));

plot(jj+1, w, 'o', 'Color', [.6 .6 .6], 'MarkerFaceColor', [.6 .6 .6]);



%%
ylabel('Weighting Index', 'FontSize', 16)
set(gca, 'FontSize', 14)
xlabel('Trial window start', 'FontSize', 16) 
set(gca, 'FontSize', 14)

ylim([-1 1])
%xlim([0 18])
set(gcf, 'color', 'w');

hh = refline(0,0);
set(hh, 'color', 'k', 'linestyle', '--');
set(gcf, 'units', 'inches', 'pos', [0 0 7 3.7])


supertitle('Nancy', 13);

wSessionsFig = gcf;
set(wSessionsFig, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 3.7]); 
figdir = '/Users/Aaron/Dropbox/twagAnalysis4.1/';
saveas(wSessionsFig, fullfile(figdir, 'wIdxChunked.pdf'))