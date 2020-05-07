subplot(233)

h=histogram(cell2mat(rval.pre(:)), 'facecolor', [.9 .9 .9]); hold on; 
h.BinWidth = 0.1;% h.Normalization = 'probability';
%errorbar(rMeans.pre.m, 400, 1.96*rMeans.pre.se, 'horizontal', 'color', [.9 .9 .9])
plot(rMeans.pre.m, 430, 'v', 'color', [.9 .9 .9])

h = histogram(cell2mat(rval.total(:)), 'facecolor', [.6 .6 0]); hold on; 
h.BinWidth = 0.1;% h.Normalization = 'probability';
%errorbar(rMeans.full.m, 400, 1.96*rMeans.full.se, 'horizontal', 'color', [.6 0 0])
plot(rMeans.full.m, 430, 'v', 'color', [.6 .6 0])

h = histogram(cell2mat(rval.post(:)), 'facecolor', [.4 .4 .4]); hold on;
h.BinWidth = 0.1; %h.Normalization = 'probability';
%errorbar(rMeans.post.m, 400, 1.96*rMeans.post.se, 'horizontal', 'color', [.4 .4 .4])
plot(rMeans.post.m, 430, 'v', 'color', [.4 .4 .4])

%legend('pre', '', 'stimulus', '', 'post', '')

title('LATE', 'fontsize', 14)

%%

subplot(234)

h = histogram(cell2mat(rval.first(:)), 'facecolor', [.5 0 0]); hold on; h.BinWidth = 0.1; h.FaceAlpha = 1;
%errorbar(rMeans.first.m, 500, 1.96*rMeans.full.se, 'horizontal', 'color', [.4 .4 0])
plot(rMeans.first.m, 250, 'v', 'color', [.5 0 0])

h = histogram(cell2mat(rval.last(:)), 'facecolor', [.75 0 0]); hold on; h.BinWidth = 0.1; h.FaceAlpha = 1;
%errorbar(rMeans.last.m, 500, 1.96*rMeans.post.se, 'horizontal', 'color', [.8 .8 0])
plot(rMeans.last.m, 250, 'v', 'color', [.75 0 0])

%legend('1st half', '', '2nd half', '')

xlabel('r_s_c', 'fontsize', 12)

%%

clr1 = clr - .2; clr1(clr1<0) = 0;
clr0 = clr - .05; clr0(clr0<0) = 0;
%clr0 = clr+.3; clr0(clr0>1) = 1;

r_full   = cell2mat(rval.total(:));
r_first   = cell2mat(rval.first(:));
r_last   = cell2mat(rval.last(:));

sameDP = cell2mat(signVal.dp(:));
sameCP = cell2mat(signVal.cp(:));
sameBoth = cell2mat(signVal.both(:));
cpPos = cell2mat(signVal.cpPos(:));
cpNeg = cell2mat(signVal.cpNeg(:));

errorbar(1, nanmean(r_full(sameCP==0)), nanstd(r_full(sameCP==0)) /sqrt(sum(sameCP==0)),...
'o', 'color', [0.5 0.5 0.5], 'lineWidth', 1.6); hold on
errorbar(2, nanmean(r_full(sameCP==1)), nanstd(r_full(sameCP==1)) /sqrt(sum(sameCP==1)),...
'o', 'color', clr, 'lineWidth', 1.6)
errorbar(3, nanmean(r_full(cpPos==1)), nanstd(r_full(cpPos==1)) /sqrt(sum(cpPos==1)),...
'o', 'color', clr1, 'lineWidth', 1.6)
errorbar(4, nanmean(r_full(cpNeg==1)), nanstd(r_full(cpNeg==1)) /sqrt(sum(cpNeg==1)),...
'o', 'color', clr0, 'lineWidth', 1.6)
xlim([0 5])
xticks(0:1:5)
xticklabels({'', 'opp', 'all same' 'same >.5', 'same<.5'})