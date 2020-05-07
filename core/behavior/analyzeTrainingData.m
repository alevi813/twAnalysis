%% analyzeTrainingData
%  get info from a twag training pds file
longTermSummary = false;

% dataPath = '/Volumes/HukLab/Macaque/Projects/leo_trainingTwag';
% pdsfile  = 'leo20191204twag.twag1346.PDS';

dataPath = '/Volumes/HukLab/Macaque/Projects/twag/leo';
pdsfile = '20191205/leo20191205twag.twag1338.PDS';

% dataPath = '/Users/aaronlevi/Desktop/leoPDS';
% pdsfile = 'leo20190806twag.twag1234.PDS';

load([dataPath filesep pdsfile], '-mat');

fprintf(['\n ' num2str(pdsfile(1:11))]);

%% stim center and theta info
[center, theta] = getDurations(PDS);

% locations(:,1) = unique(center(:,1));
% locations(:,2) = unique(center(:,2));

angles = unique(theta);

%% n goodtrials / out of all

goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
goodtrial = cell2mat(goodtrial);

stimDistNum = cellfun(@(x) x.theseGabs.stimDistNum, PDS.data, 'UniformOutput', false);
stimDistNum = cell2mat(stimDistNum);
stimDistNum = stimDistNum(goodtrial==1);

sumCoh = cellfun(@(x) x.theseGabs.sumCoh, PDS.data, 'UniformOutput', false);
sumCoh = cell2mat(sumCoh);
sumCoh = sumCoh(goodtrial==1);

for ii = 1:length(PDS.data)
    if ~isfield(PDS.data{ii}, 'targRchosen')
        PDS.data{ii}.targRchosen = nan;
    end
    if ~isfield(PDS.data{ii}, 'correct')
        PDS.data{ii}.correct = nan;
    end
    if ~isfield(PDS.data{ii}.theseGabs, 'targCorrect')
        PDS.data{ii}.theseGabs.targCorrect = nan;
    end
end

%fprintf(['\n center = ' num2str(locations)])
fprintf(['\n theta = ' num2str(angles')])
fprintf(['\n n goodtrials = ' num2str(nansum(goodtrial)) ' of ' num2str(length(PDS.data)) ' trials']);

%% n two target trials / all good
% everything below will only be calculated on two targ trials.

% if isfield(PDS.data{1}.stimulus, 'oneTarg')
%     oneTarg     = cell2mat(cellfun(@(x) x.stimulus.oneTarg, PDS.data, 'UniformOutput', false));
%     oneTarg     = oneTarg(goodtrial==1);
% else 
%     oneTarg = NaN;
% end

% fprintf(['\n n two target trials = ' num2str(nansum(~oneTarg)) ' of ' num2str(nansum(goodtrial)) ' trials']);


%% prop R correct

targRcorrect = cellfun(@(x) x.targRcorrect, PDS.data, 'UniformOutput', false);
targRcorrect = targRcorrect(goodtrial==1);
%targRcorrect = targRcorrect(~oneTarg);
targRcorrect = cell2mat(targRcorrect);

propTargRcorrect   = sum(targRcorrect)/length(targRcorrect);

fprintf(['\n prop. R correct = ' num2str(propTargRcorrect)]);
%% prop R chosen
targRchosen = cellfun(@(x) x.targRchosen, PDS.data, 'UniformOutput', false);
targRchosen = targRchosen(goodtrial==1);
%targRchosen = targRchosen(~oneTarg);
targRchosen = cell2mat(targRchosen);

propTargRchosen   = sum(targRchosen)/length(targRchosen);

fprintf(['\n prop. R choices = ' num2str(propTargRchosen)]);

%% total pc

correct  = cellfun(@(x) x.correct, PDS.data, 'UniformOutput', false);
correct  = cell2mat(correct(goodtrial==1));
%correct  = correct(~oneTarg);

pc = sum(correct)/length(correct);

fprintf(['\n total pc = ' num2str(pc)])

% get sumcoh as well
%sumCoh = sumCoh(~oneTarg);

%% proportion global motion

% global motion trials
for itrial = 1:length(PDS.data)
    gmExist(itrial) = isfield(PDS.data{itrial}.stimulus, 'globalMotion');
    
    if gmExist(itrial) & PDS.data{itrial}.stimulus.globalMotion.trialFlag == 1
        gm.flag(itrial) = 1;
%        gm.mag(itrial)  = PDS.data{itrial}.stimulus.globalMotion.magnitude;
    else
        gm.flag(itrial) = 0;
%        gm.mag(itrial)  = 0;
    end
end
gm.flag = gm.flag(goodtrial==1);
%gm.flag = gm.flag(~oneTarg);

%propgm = sum(gm.flag)/length(gm.flag);

fprintf(['\n\n n global motion = ' num2str(sum(gm.flag)) ' of ' num2str(length(gm.flag)) ]);


%% global motion pc
% gm percent correct
gm.correct = correct(gm.flag==1);
gm.pc = sum(gm.correct)/length(gm.correct);

%fprintf(['\n all global motion pc = ' num2str(gm.pc) ' (' num2str(sum(gm.flag)) ' trials)'  ])

% global motion pc, by stimDist
%gm.stimDistNum = stimDistNum(~oneTarg);
gm.stimDistNum = stimDistNum;
gm.stimDistNum = gm.stimDistNum(gm.flag==1);

gm.distHi  = gm.stimDistNum==1 | gm.stimDistNum==5;
gm.distMed = gm.stimDistNum==2 | gm.stimDistNum==4;
gm.distLo  = gm.stimDistNum==3;

gm.pcHi  = sum(gm.correct(gm.distHi))   / length(gm.correct(gm.distHi));
gm.pcMed = sum(gm.correct(gm.distMed)) / length(gm.correct(gm.distMed));
gm.pcLo  = sum(gm.correct(gm.distLo)) / length(gm.correct(gm.distLo));

%fprintf(['\n strong global motion pc = ' num2str(gm.pcHi) ' (' num2str(sum(gm.distHi)) ' trials) '])
%fprintf(['\n weak global motion pc = ' num2str(gm.pcMed) ' (' num2str(sum(gm.distMed)) ' trials)'])
%fprintf(['\n zero-mean global motion pc = ' num2str(gm.pcLo) ' (' num2str(sum(gm.distLo)) ' trials) \n'])

%% local motion pc

lm.correct = correct(gm.flag==0);
lm.pc = sum(lm.correct)/length(lm.correct);

fprintf(['\n n local motion = ' num2str(sum(~gm.flag)) ' of ' num2str(length(~gm.flag)) ]);
fprintf(['\n all local motion pc = ' num2str(lm.pc) ' (' num2str(sum(~gm.flag)) ' trials)'  ])

% global motion pc, by stimDist
%lm.stimDistNum = stimDistNum(~oneTarg);
lm.stimDistNum = stimDistNum;
lm.stimDistNum = lm.stimDistNum(~gm.flag);

lm.distHi  = lm.stimDistNum==1 | lm.stimDistNum==5;
lm.distMed = lm.stimDistNum==2 | lm.stimDistNum==4;
lm.distLo  = lm.stimDistNum==3;

lm.pcHi  = sum(lm.correct(lm.distHi))   / length(lm.correct(lm.distHi));
lm.pcMed = sum(lm.correct(lm.distMed)) / length(lm.correct(lm.distMed));
lm.pcLo  = sum(lm.correct(lm.distLo)) / length(lm.correct(lm.distLo));

isFro  = stimDistNum == 6;
nFro   = sum(isFro);
froCho = targRchosen(isFro);
propFroCho = sum(froCho) / nFro;

fprintf(['\n strong local motion pc = ' num2str(lm.pcHi) ' (' num2str(sum(lm.distHi)) ' trials) '])
fprintf(['\n weak local motion pc = ' num2str(lm.pcMed) ' (' num2str(sum(lm.distMed)) ' trials)'])
fprintf(['\n zero-mean local motion pc = ' num2str(lm.pcLo) ' (' num2str(sum(lm.distLo)) ' trials) \n'])

fprintf(['\n frozen trials prop R choices = ' num2str(propFroCho) ' (' num2str(nFro) ' trials) \n'])

%% plot
%  running average

% n trial sliding window
window = 20; 

for ii = 1:length(correct)-window

    tempCorrect     = correct(ii:ii+window);
    propCorrect(ii) = sum(tempCorrect) / length(tempCorrect);
    
    tempChoice     = targRchosen(ii:ii+window);
    propChoice(ii) = sum(tempChoice) / length(tempChoice);
    
    tempTargRcorrect     = targRcorrect(ii:ii+window);
    propTargRcorrect(ii) = sum(tempTargRcorrect) / length(tempTargRcorrect);

    tempSumCoh = abs(sumCoh(ii:ii+window));
    avgSumCoh(ii)  = (sum(tempSumCoh) / length(tempSumCoh)) / 100;
end


figure
subplot(3,1,1); hold on
title('single trial')
plot(abs(sumCoh));
ylabel('coherence')

subplot(3,1,2); hold on
title('sliding trial window = 20')
plot([0 length(propCorrect)], [0.5 0.5], 'k--')
plot(propCorrect);
plot(avgSumCoh, 'g--')
legend('', 'prop. correct', 'prop. coh')

subplot(3,1,3); hold on
plot([0 length(propChoice)], [0.5 0.5], 'k--')
plot(propChoice);
plot(propTargRcorrect);
legend('', 'prop. r choices', 'prop. r correct')
xlabel('trial')

supertitle(num2str(pdsfile(1:11)), 13);

% imagesc the choice and the correct answer
figure
imagesc( [targRchosen ; targRcorrect]' ); colorbar
set(gca, 'XTick', [1 2], 'XTickLabel', {'Choice', 'Correct'}, 'FontSize', 14)


%% pmf + ppk
trStart = 1;
pp = analyzeBehaviorSession(dataPath, pdsfile, trStart, 1);
getThreshold(pp.pmfUnfolded, .75);
%figure(4); 
subplot(121);
text( -50, .75, ['R thresh = ' num2str(pp.pmfUnfolded.threshValue)])

getThreshold(pp.pmfUnfolded, .25);
text( 5, .25, ['L thresh = ' num2str(pp.pmfUnfolded.threshValue)])




%%
if longTermSummary
% training summary
dataPath = '/Volumes/LKCLAB/EPHYS/DATA/Projects/leo_trainingTwag';

files = dir(dataPath);
files = files(end-19:end);

for iFile = 1:length(files)
load([dataPath filesep files(iFile).name], '-mat');

% ntrials
goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
goodtrial = cell2mat(goodtrial);
nGood(iFile) = nansum(goodtrial);

% center & theta -- get from recreateParams
[center, theta] = getDurations(PDS);

locations{iFile}(:,1) = unique(center(:,1));
locations{iFile}(:,2) = unique(center(:,2));

angles{iFile} = unique(theta);

% pmf threshold
for ii = 1:length(PDS.data)
    if ~isfield(PDS.data{ii}, 'targRchosen')
        PDS.data{ii}.targRchosen = nan;
    end
    if ~isfield(PDS.data{ii}, 'correct')
        PDS.data{ii}.correct = nan;
    end
    if ~isfield(PDS.data{ii}.theseGabs, 'targCorrect')
        PDS.data{ii}.theseGabs.targCorrect = nan;
    end
end

sumCoh = cellfun(@(x) x.theseGabs.sumCoh, PDS.data, 'UniformOutput', false);
sumCoh = cell2mat(sumCoh);
sumCoh = sumCoh(goodtrial==1);

targRchosen = cellfun(@(x) x.targRchosen, PDS.data, 'UniformOutput', false);
targRchosen = targRchosen(goodtrial==1);
targRchosen = cell2mat(targRchosen);

pmf   = pmfTools(sumCoh, targRchosen, 'nBins', 12);
pmf   = fit(pmf);

pmf   = getThreshold(pmf, .75);
rThresh = pmf.threshValue;
pmf   = getThreshold(pmf, .25);
lThresh = pmf.threshValue;

thresh(iFile) = nanmean([abs(rThresh) abs(lThresh)]);



end % file loop

end