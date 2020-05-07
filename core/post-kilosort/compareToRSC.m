%% Noise corr, cp, and d'

% set up some paths
comp = getComp;

if strcmp(comp, 'laptop')
    dataPath{1} = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/nancy';
    dataPath{2} = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo';
else
    dataPath{1} = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy';
    dataPath{2} = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo';
end

%condition = {'flat', 'late', 'early'}; % flat/late/early --- use 'good' for all sessions
condition = {'late'}; % flat/late/early --- one condition at a time.

% load neurons
S = loadNeuronsByCondition(condition);

for iCond = 1:length(condition)
    % set up some colors
    switch condition{iCond}
        case 'flat'
            clr = [19/255 40/255 230/255];
            clr = [0 0 .6];
        case 'late'
            clr = [255/255 234/255 0/255];
            clr = [.6 .6 0];
        case 'early'
            clr = [214/255 4/255 0/255];
            clr = [.6 0 0];
    end
    
    % load dprime
    allVal.dp = arrayfun(@(x) x.model(1).dprime, S);

    % get cp (full stimulus window)
    allVal.cp = grandCP_test(S, dataPath);
    
    % CORRELATIONS
    doShuffle = false;
    [rval, signVal, pairDiff, pairMean] = noiselCorrelations(S, allVal, dataPath, doShuffle);

end

%%
sameDP = logical(cell2mat(signVal.dp(:)));

tmp = cell2mat(rval.total(:));
tmpFirst = cell2mat(rval.first(:));
tmpLast  = cell2mat(rval.last(:));

rMeans.sameDP.full.m = nanmean(tmp(sameDP));
rMeans.sameDP.full.se  = nanstd(tmp(sameDP)) / sqrt(length(tmp(sameDP)));

rMeans.sameDP.first.m = nanmean(tmpFirst(sameDP));
rMeans.sameDP.first.se  = nanstd(tmpFirst(sameDP)) / sqrt(length(tmpFirst(sameDP)));

rMeans.sameDP.last.m = nanmean(tmpLast(sameDP));
rMeans.sameDP.last.se  = nanstd(tmpLast(sameDP)) / sqrt(length(tmpLast(sameDP)));

rMeans.full.m   = nanmean(cell2mat(rval.total(:))); 
rMeans.full.se  = nanstd(cell2mat(rval.total(:))) / sqrt(length(cell2mat(rval.total(:))));

rMeans.first.m   = nanmean(cell2mat(rval.first(:)));
rMeans.first.se  = nanstd(cell2mat(rval.first(:))) / sqrt(length(cell2mat(rval.first(:))));

rMeans.last.m     = nanmean(cell2mat(rval.last(:)));
rMeans.last.se  = nanstd(cell2mat(rval.last(:))) / sqrt(length(cell2mat(rval.last(:))));

rMeans.pre.m   = nanmean(cell2mat(rval.pre(:)));
rMeans.pre.se  = nanstd(cell2mat(rval.pre(:))) / sqrt(length(cell2mat(rval.pre(:))));

rMeans.post.m   = nanmean(cell2mat(rval.post(:)));
rMeans.post.se  = nanstd(cell2mat(rval.pre(:))) / sqrt(length(cell2mat(rval.pre(:))));


%% plot

% pre, during, post stimulus rsc for all conditions
figure(1); 
% subplot(121); hold on
% errorbar([1 2 3], [rMeans.pre.m rMeans.full.m rMeans.post.m],...
%     [rMeans.pre.se rMeans.full.se rMeans.post.se], '-o', 'color', clr)
% axis square
% xlim([0 4])
% xticklabels({'', 'pre', 'stimulus', 'post', ''})

subplot(236); hold on
% errorbar([1 2], [rMeans.first.m rMeans.last.m],...
%     [rMeans.first.se  rMeans.last.se], '-o', 'color', [.6 .6 .6])

errorbar([1 2], [rMeans.sameDP.first.m rMeans.sameDP.last.m],...
    [rMeans.sameDP.first.se  rMeans.sameDP.last.se], '-o', 'color', [.5 .5 0])
axis square
xlim([0 3])
xticklabels({'', 'first', 'last', ''})
% errorbar([1 2 3 4], [rMeans.pre.m rMeans.first.m rMeans.last.m rMeans.post.m],...
%     [rMeans.pre.se rMeans.first.se  rMeans.last.se rMeans.post.se], '-o', 'color', clr)
% axis square
% xlim([0 5])
% xticklabels({'', 'pre', 'first', 'last', 'post', ''})

%%

subplot(232)
%plot(rMeans.pre.m, 430, 'v', 'color', clr)

h = histogram(cell2mat(rval.total(:)), 'facecolor', [.6 .6 .6]); hold on; 
h.BinWidth = 0.1;% h.Normalization = 'probability';
%errorbar(rMeans.full.m, 400, 1.96*rMeans.full.se, 'horizontal', 'color', [.6 0 0])
%plot(rMeans.full.m, 430, 'v', 'color', clr)


% %%
% % rsc by same/different sign of dp cp
% figure(2)
% r_full   = cell2mat(rval.total(:));
% r_first   = cell2mat(rval.first(:));
% r_last   = cell2mat(rval.last(:));
% 
% sameDP = cell2mat(signVal.dp(:));
% sameCP = cell2mat(signVal.cp(:));
% sameBoth = cell2mat(signVal.both(:));
% 
% subplot(131); hold on
% clr1 = clr - .2; clr1(clr1<0) = 0;
% %clr0 = clr - .05; clr0(clr0<0) = 0;
% clr0 = clr+.3; clr0(clr0>1) = 1;
% errorbar(1, nanmean(r_full(sameDP==1)), nanstd(r_full(sameDP==1)) /sqrt(sum(sameDP==1)),...
%     'o', 'color', clr1, 'lineWidth', 1.6)
% hold on
% errorbar(1.5, nanmean(r_full(sameDP==0)), nanstd(r_full(sameDP==0)) /sqrt(sum(sameDP==0)),...
%     'o', 'color', clr0, 'lineWidth', 1.6)
% 
% subplot(132); hold on
% errorbar(1, nanmean(r_full(sameCP==1)), nanstd(r_full(sameCP==1)) /sqrt(sum(sameCP==1)),...
%     'o', 'color', clr1, 'lineWidth', 1.6)
% hold on
% errorbar(1.5, nanmean(r_full(sameCP==0)), nanstd(r_full(sameCP==0)) /sqrt(sum(sameCP==0)),...
%     'o', 'color', clr0, 'lineWidth', 1.6)
% 
% subplot(133); hold on
% errorbar(1, nanmean(r_full(sameBoth==1)), nanstd(r_full(sameBoth==1)) /sqrt(sum(sameBoth==1)),...
%     'o', 'color', clr1, 'lineWidth', 1.6)
% hold on
% errorbar(1.5, nanmean(r_full(sameBoth==0)), nanstd(r_full(sameBoth==0)) /sqrt(sum(sameBoth==0)),...
%     'o', 'color', clr0, 'lineWidth', 1.6)
% 
% % xlim([-1 3.5])
% % xticks(-1:.5:3.5)
% % xticklabels({'', 'same', 'opp', '', 'same', 'opp', '', 'same', 'opp', ''})
% 
% %% scatterplots
% 
% nixr = isnan(r_full);
% clr = clr-.2; clr(clr<0) = 0;
% 
% figure(3)
% subplot(232)
% dpDiff = cell2mat(pairDiff.dp(:));
% plot(abs(dpDiff), r_full, 'o', 'color', clr)
% % plot(abs(dpDiff), r_first, 'ko'); hold on
% % plot(abs(dpDiff), r_last, '^', 'color', [.6 .6 .6])
% 
% xlabel('dprime difference')
% ylabel('r_s_c')
% nix = isnan(dpDiff);
% tmpR = corrcoef(dpDiff(~nix&~nixr), r_full( (~nix&~nixr)));
% title(['r = ' num2str(tmpR(1,2))])
% 
% subplot(235)
% dpMean = cell2mat(pairMean.dp(:));
% plot(dpMean, r_full, 'o', 'color', clr)
% % plot(dpMean, r_first, 'ko'); hold on
% % plot(dpMean, r_last, '^', 'color', [.6 .6 .6])
% 
% xlabel('dprime mean')
% ylabel('r_s_c')
% nix = isnan(dpMean);
% tmpR = corrcoef(dpMean(~nix&~nixr), r_full( (~nix&~nixr)));
% title(['r = ' num2str(tmpR(1,2))])
% 
% figure(4)
% subplot(233)
% cpDiff = cell2mat(pairDiff.cp(:));
% plot(abs(cpDiff), r_full, 'o', 'color', clr)
% % plot(abs(cpDiff),  r_first, 'ko'); hold on
% % plot(abs(cpDiff), r_last, '^', 'color', [.6 .6 .6])
% 
% xlabel('cp difference')
% ylabel('r_s_c')
% nix = isnan(cpDiff);
% tmpR = corrcoef(cpDiff(~nix&~nixr), r_full( (~nix&~nixr)));
% title(['r = ' num2str(tmpR(1,2))])
% 
% subplot(236)
% cpMean = cell2mat(pairMean.cp(:));
% plot(cpMean, r_full, 'o', 'color', clr)
% % plot(cpMean, r_first, 'ko'); hold on
% % plot(cpMean, r_last, '^', 'color', [.6 .6 .6])
% 
% xlabel('cp mean')
% ylabel('r_s_c')
% nix = isnan(cpMean);
% tmpR = corrcoef(cpMean(~nix&~nixr), r_full( (~nix&~nixr)));
% title(['r = ' num2str(tmpR(1,2))])
% 
% 
