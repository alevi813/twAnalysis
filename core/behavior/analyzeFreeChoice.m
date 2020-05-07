% alvin free choice training
% 
% date =  '01/11/2018';
% pdsAll{1} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180111twag.twag1233.PDS';
% pdsAll{2} ='/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180111twag.twag1257.PDS';
% pdsAll{3} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180111twag.twag1326.PDS';
% pdsAll{4} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180111twag.twag1421.PDS';

% date =  '01/18/2018';
% pdsAll{1} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180118twag.twag1151.PDS';
% pdsAll{2} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180118twag.twag1207.PDS';
% pdsAll{3} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180118twag.twag1222.PDS';
% pdsAll{4} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180118twag.twag1343.PDS';
% pdsAll{5} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180118twag.twag1400.PDS';

%date =  '01/19/2018';
%pdsAll{1} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180119twag.twag1208.PDS';
 
%pdsAll{1} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/gaborMotion/alvin20180119twag.twag1424.PDS';

date =  '01/23/2018';
pdsAll{1} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180123twag.twag1127.PDS';
pdsAll{2} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/freeChoice/alvin20180123twag.twag1158.PDS';

%pdsAll{1} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/gaborMotion/alvin20180123twag.twag1214.PDS';
%pdsAll{1} = '/Volumes/LKCLAB/EPHYS/DATA/Projects/alvin_trainingTwag/gaborMotion/alvin20180123twag.twag1308.PDS';

choiceAll = [];
rewardprob = nan(length(pdsAll), 1);

for ii = 1:length(pdsAll)
    
    load(pdsAll{ii}, '-mat');
    if isfield(PDS.initialParametersMerged.behavior.reward, 'prob')
       rewardprob(ii) = PDS.initialParametersMerged.behavior.reward.prob;
    end
    
    for ii = 1:length(PDS.data)
        if ~isfield(PDS.data{ii}, 'targRchosen')
            PDS.data{ii}.targRchosen = nan;
        end
    end
    
    goodtrial = cell2mat(cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false));
    goodtrial(isnan(goodtrial)) = 0;
    
    targRchosen = cell2mat(cellfun(@(x) x.targRchosen, PDS.data, 'UniformOutput', false));
    targRchosen = targRchosen(goodtrial==1);
    
    choiceAll = [choiceAll, targRchosen];
end

propR = sum(choiceAll)/length(choiceAll);
propL = 1-propR;

figure;
bar([propL, propR])
title(date);

ylim([0 1])
xlabel('Choice')
ylabel('proportion')
set(gca, 'XTickLabel', {'left'; 'right'} );

% histogram(choiceAll);