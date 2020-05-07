function [allrevco, frozenonly] = unitCP_postKilo(centerTimeGT, pp, theseSpikes, dirPref)

% calculate cp time course
% inputs: 
%     theseSpikes  - spike times for specific unit
%     pp           - post-processed version of PDS
%     centerTimeGT - centering time for goodtrials (ie motion on)
%     dirPref      - direction preference of unit - 1 right, 0 left
% outputs:
%     cp struct with fields for allrevco and frozenonly. Each of which was
%     mean and std fields (nUnits x nTimeBins)


%% allrevco
if dirPref == 1
    T1ix=pp.targRchosen(:)==1 & pp.stimDistNum(:)==3; %t1ix if targ 1 (right is preferred dir)
    T2ix=pp.targRchosen(:)==0 & pp.stimDistNum(:)==3;
else
    T2ix=pp.targRchosen(:)==1 & pp.stimDistNum(:)==3;
    T1ix=pp.targRchosen(:)==0 & pp.stimDistNum(:)==3;
end

if length(T1ix) ~= length(centerTimeGT)   % fix for sessions that run longer after turning of recording. Must vet this.
    T1ix = T1ix(1:length(centerTimeGT));
    T2ix = T2ix(1:length(centerTimeGT));
end


% align
[mT1,sT1,bc, ~, trSpkCntT1]=pdsa.eventPsth(theseSpikes, centerTimeGT(T1ix), [-.5 1.5], .001, ones(100,1)/100);
[mT2,sT2, ~, ~, trSpkCntT2]=pdsa.eventPsth(theseSpikes, centerTimeGT(T2ix), [-.5 1.5], .001, ones(100,1)/100);

allTrSpikes = [trSpkCntT1; trSpkCntT2];
allTrChoice = [ones(size(trSpkCntT1, 1),1); zeros(size(trSpkCntT2, 1),1)];

allrevco.nTrials(1,1) = size(allTrSpikes, 1);

% ajl debugging -- removing trials with 0 or nan spk counts
nospikes = nan(size(allTrSpikes));
for itrial = 1:size(allTrSpikes, 1)
    if allTrSpikes(itrial, :) == 0
        nospikes(itrial,:) = 1;
    else
        nospikes(itrial,:) = 0;
    end
end
%badTrials    = allTrSpikes == 0 | isnan(allTrSpikes);
badTrials    = nospikes | isnan(allTrSpikes);
allTrSpikes  = allTrSpikes(~badTrials(:,1),:);
allTrChoice  = allTrChoice(~badTrials(:,1));

allrevco.nTrials(1,2) = size(allTrSpikes, 1);

% calculate cp
[allrevco.m, allrevco.s] = choiceProbabilityCalculate(allTrSpikes, allTrChoice);

allrevco.grandcp = mean(allrevco.m(600:1600)); % HARDCODED, BE CAREFUL. PROBALY WANT TO USE FUNCTION 'GRANDCP' INSTEAD
% disp(['went from ' num2str(tmp_ntrials) ' to ' num2str(allrevco.nTrials) ' revco trials' ]);

%% frozen only
if dirPref == 1
    T1ix=pp.targRchosen(:)==1 & pp.stimDistNum(:)==6; %t1ix if targ 1 (right is preferred dir)
    T2ix=pp.targRchosen(:)==0 & pp.stimDistNum(:)==6;
else
    T2ix=pp.targRchosen(:)==1 & pp.stimDistNum(:)==6;
    T1ix=pp.targRchosen(:)==0 & pp.stimDistNum(:)==6;
end

if length(T1ix) ~= length(centerTimeGT)   % fix for sessions that run longer after turning of recording. Must vet this.
    T1ix = T1ix(1:length(centerTimeGT));
    T2ix = T2ix(1:length(centerTimeGT));
end

% align
[mT1,sT1,bc, ~, trSpkCntT1]=pdsa.eventPsth(theseSpikes, centerTimeGT(T1ix), [-.5 1.5], .001, ones(100,1)/100);
[mT2,sT2, ~, ~, trSpkCntT2]=pdsa.eventPsth(theseSpikes, centerTimeGT(T2ix), [-.5 1.5], .001, ones(100,1)/100);

allTrSpikes = [trSpkCntT1; trSpkCntT2];
allTrChoice = [ones(size(trSpkCntT1, 1),1); zeros(size(trSpkCntT2, 1),1)];

frozenonly.nTrials(1,1) = size(allTrSpikes, 1);

% ajl debugging -- removing trials with 0 or nan spk counts
nospikes = nan(size(allTrSpikes));
for itrial = 1:size(allTrSpikes, 1)
    if allTrSpikes(itrial, :) == 0
        nospikes(itrial,:) = 1;
    else
        nospikes(itrial,:) = 0;
    end
end
%badTrials    = allTrSpikes == 0 | isnan(allTrSpikes);
badTrials    = nospikes | isnan(allTrSpikes);
allTrSpikes  = allTrSpikes(~badTrials(:,1),:);
allTrChoice  = allTrChoice(~badTrials(:,1));

frozenonly.nTrials(1,2) = size(allTrSpikes, 1);

% calculate cp
[frozenonly.m, frozenonly.s] = choiceProbabilityCalculate(allTrSpikes, allTrChoice);

% frozenonly.grandcp = mean(frozenonly.m(600:1600)); % HARDCODED, BE CAREFUL. PROBALY WANT TO USE FUNCTION 'GRANDCP' INSTEAD
% disp(['went from ' num2str(tmp_ntrials) ' to ' num2str(frozenonly.nTrials) ' frozen trials']);
end