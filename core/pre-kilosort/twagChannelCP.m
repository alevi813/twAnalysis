function [cp, bc, allTrSpikes, allTrChoice] = twagChannelCP(cpPlot, channelNumber, spikeTimes, pp, centerTime, goodtrial, dirPref)

% function to calculate and plot the time course of choice probability for
% twag experiments
%
% INPUTS:
%    - channelNumber: vector of channel numbers from recording that you
%    wish to load.
%    - spikeTimes: nCh x 1 cell with spike times from plx_ts
%    - pp: the processed behavior file
%    - centerTime: beginning of window for PSTH using trial start and
%    motion onset times, calculated in twagPlxAnalysis
%    - goodtrial: index of completed 'goodtrials' from PDS, renamed in
%    twagPlxAnalysis
%    - dirPref: direction prefence (string), from twagPSTH
% OUTPUTS:
%    - grandCP_allRevco: overall choice probability calculated from
%    all revco trials.
%    - grandCP_frozen: overall of choice probability calculated
%    from frozen trials only.
% -------------------------------------------------------------------------


nCh         = length(channelNumber);
allTrSpikes = cell(nCh, 1);
allTrChoice = cell(nCh, 1);

mCPallRevco      = cell(nCh, 1);
sCPallRevco      = cell(nCh, 1);
grandCP_allRevco = nan(nCh, 1);

mCPfrozen      = cell(nCh, 1);
sCPfrozen      = cell(nCh, 1);
grandCP_frozen = nan(nCh, 1);

%% Calculate time course of CP using all revco trials
if cpPlot == 1
    figure; hold on
end

kk = 1;
%for jj = 1:nCh
for jj = channelNumber(kk):channelNumber(nCh)
    
    if length(goodtrial)==length(centerTime)     % fix for sessions that run longer after turning of recording. Must vet this.
        centerTimeGT = centerTime(goodtrial==1);
    else
        goodtrial = goodtrial(1:length(centerTime));
        centerTimeGT = centerTime(goodtrial==1);
    end
    
    % assign indices according to directional preference
    if dirPref(jj) == 1
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
    
    [m1,s1,bc, ~, trSpkCnt1]=pdsa.eventPsth(spikeTimes{jj,1}, centerTimeGT(T1ix), [-.2 1.5], .001, ones(100,1)/100);
    [m2,s2, ~, ~, trSpkCnt2]=pdsa.eventPsth(spikeTimes{jj,1}, centerTimeGT(T2ix), [-.2 1.5], .001, ones(100,1)/100);
    
    %         figure
    %         plot(bc, m1, 'b', bc, m2, 'r', ...
    %             bc, m1+s1, 'b--', bc, m1-s1, 'b--', ...
    %             bc, m2+s2, 'r--', bc, m2-s2, 'r--')
    %     legend('Choice R', 'Choice L');
    
    %time course of CP
    allTrSpikes{jj,1} = [trSpkCnt1; trSpkCnt2];
    allTrChoice{jj,1} = [ones(size(trSpkCnt1, 1),1); zeros(size(trSpkCnt2, 1),1)];
    %allTrChoice{jj,1} = [ones(size(trSpkCnt2, 1),1); zeros(size(trSpkCnt1, 1),1)];
    
    [mCPallRevco{jj,1}, sCPallRevco{jj,1}] = choiceProbabilityCalculate(allTrSpikes{jj,1}, allTrChoice{jj,1});
    
    if cpPlot == 1
        %subplot(1, nCh, kk); hold on
        %subplot(6, 4, kk); hold on
        subplot(4, 8, kk); hold on
        %refline(0, .5);
        plot(bc, mCPallRevco{jj,1},'k',...
            bc, mCPallRevco{jj,1}+sCPallRevco{jj,1}, 'k--',...
            bc, mCPallRevco{jj,1}-sCPallRevco{jj,1}, 'k--');
        xlim([0 1.5]);
    end
    
    %grandCP_allRevco(jj) = mean(mCPallRevco{jj,1}(350:1300));
    grandCP_allRevco(jj) = mean(mCPallRevco{jj,1}(350:1400));
    
    if cpPlot == 1
        text(0.1, max(mCPallRevco{jj,1}), num2str(round(grandCP_allRevco(jj),2)));
    end
    
    kk = kk+1;
end

if cpPlot == 1
    supertitle('CP from all RevCo trials', 12);
end
%% Calculate time course of CP using only frozen trials
if cpPlot == 1
    figure; hold on
end

kk = 1;
%for jj = 1:nCh
for jj = channelNumber(kk):channelNumber(nCh)
    
    % assign directional indicees
    if dirPref(jj) == 1
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
    
    [m1,s1,bc, ~, trSpkCnt1]=pdsa.eventPsth(spikeTimes{jj,1}, centerTimeGT(T1ix), [-.2 1.5], .001, ones(100,1)/100);
    [m2,s2, ~, ~, trSpkCnt2]=pdsa.eventPsth(spikeTimes{jj,1}, centerTimeGT(T2ix), [-.2 1.5], .001, ones(100,1)/100);
    
    %     figure
    %     plot(bc, m1, 'b', bc, m2, 'r', ...
    %         bc, m1+s1, 'b--', bc, m1-s1, 'b--', ...
    %         bc, m2+s2, 'r--', bc, m2-s2, 'r--')
    %     legend('Choice R', 'Choice L');
    
    %time course of CP
    allTrSpikes{jj,1} = [trSpkCnt1; trSpkCnt2];
    allTrChoice{jj,1} = [ones(size(trSpkCnt1, 1),1); zeros(size(trSpkCnt2, 1),1)];
    %allTrChoice{jj,1} = [ones(size(trSpkCnt2, 1),1); zeros(size(trSpkCnt1, 1),1)];
    
    [mCPfrozen{jj,1}, sCPfrozen{jj,1}] = choiceProbabilityCalculate(allTrSpikes{jj,1}, allTrChoice{jj,1});
    
    if cpPlot == 1
        %subplot(1, nCh, kk); hold on
        subplot(6, 4, kk); hold on
        %refline(0, .5);
        plot(bc, mCPfrozen{jj,1},'k',...
            bc, mCPfrozen{jj,1}+sCPfrozen{jj,1}, 'k--',...
            bc, mCPfrozen{jj,1}-sCPfrozen{jj,1}, 'k--');
        xlim([0 1.5]);
    end
    
    %grandCP_frozen(jj) = mean(mCPfrozen{jj,1}(350:1300));
    grandCP_frozen(jj) = mean(mCPfrozen{jj,1}(350:1400));
    
    if cpPlot == 1
        text(0.1, max(mCPfrozen{jj,1}), num2str(round(grandCP_frozen(jj),2)));
    end
    
    kk = kk+1;
end

if cpPlot == 1
    supertitle('CP from frozen trials only', 12);
end

%% make a nice cp struct 
cp.channel.allRevco.m     = mCPallRevco;
cp.channel.allRevco.s     = sCPallRevco;
cp.channel.allRevco.grand = grandCP_allRevco;

cp.channel.frozen.m     = mCPfrozen;
cp.channel.frozen.s     = sCPfrozen;
cp.channel.frozen.grand = grandCP_frozen;

%%
% %% Calculate time course of CP using subthreshold
% figure; hold on
% kk = 1;
% %for jj = 1:nCh
% for jj = channelNumber(kk):channelNumber(nCh)
%
%     pp.pmfUnfolded.getThreshold(.75);
%
%     % assign directional indicees
%     if dirPref(jj) == 1
%         T1ix=pp.targRchosen(:)==1 & pp.stimDistNum(:)==3 & (abs(pp.sumCoh) < pp.pmfUnfolded.threshValue)'; %t1ix if targ 1 (right is preferred dir)
%         T2ix=pp.targRchosen(:)==0 & pp.stimDistNum(:)==3 & (abs(pp.sumCoh) < pp.pmfUnfolded.threshValue)';
%     else
%         T2ix=pp.targRchosen(:)==1 & pp.stimDistNum(:)==3 & (abs(pp.sumCoh) < pp.pmfUnfolded.threshValue)';
%         T1ix=pp.targRchosen(:)==0 & pp.stimDistNum(:)==3 & (abs(pp.sumCoh) < pp.pmfUnfolded.threshValue)';
%     end
%
%     if length(T1ix) ~= length(centerTimeGT)   % fix for sessions that run longer after turning of recording. Must vet this.
%         T1ix = T1ix(1:length(centerTimeGT));
%         T2ix = T2ix(1:length(centerTimeGT));
%     end
%
%     [m1,s1,bc, ~, trSpkCnt1]=pdsa.eventPsth(spikeTimes{jj,1}, centerTimeGT(T1ix), [-.2 1.5], .001, ones(100,1)/100);
%     [m2,s2, ~, ~, trSpkCnt2]=pdsa.eventPsth(spikeTimes{jj,1}, centerTimeGT(T2ix), [-.2 1.5], .001, ones(100,1)/100);
%
%     %     figure
%     %     plot(bc, m1, 'b', bc, m2, 'r', ...
%     %         bc, m1+s1, 'b--', bc, m1-s1, 'b--', ...
%     %         bc, m2+s2, 'r--', bc, m2-s2, 'r--')
%     %     legend('Choice R', 'Choice L');
%
%     %time course of CP
%     allTrSpikes{jj,1} = [trSpkCnt1; trSpkCnt2];
%     allTrChoice{jj,1} = [ones(size(trSpkCnt1, 1),1); zeros(size(trSpkCnt2, 1),1)];
%     %allTrChoice{jj,1} = [ones(size(trSpkCnt2, 1),1); zeros(size(trSpkCnt1, 1),1)];
%
%     [mCPsomeRevco{jj,1}, sCPsomeRevco{jj,1}] = choiceProbabilityCalculate(allTrSpikes{jj,1}, allTrChoice{jj,1});
%
%     %figure; hold on
%     %subplot(2, nCh, jj+nCh); hold on
%     subplot(1, nCh, kk); hold on
%     %refline(0, .5);
%     plot(bc, mCPsomeRevco{jj,1},'k',...
%         bc, mCPsomeRevco{jj,1}+sCPsomeRevco{jj,1}, 'k--',...
%         bc, mCPsomeRevco{jj,1}-sCPsomeRevco{jj,1}, 'k--');
%     xlim([0 1.5]);
%
% %    grandCP_someRevco(jj) = mean(mCPsomeRevco{jj,1}(350:1300));
%     grandCP_someRevco(jj) = mean(mCPsomeRevco{jj,1}(350:1400));
%     text(0.1, max(mCPsomeRevco{jj,1}), num2str(round(grandCP_someRevco(jj),2)));
%
%     kk = kk+1;
% end
% supertitle('CP from subthreshold RevCo trials', 12);
