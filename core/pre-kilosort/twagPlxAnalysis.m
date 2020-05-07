%% script to calculate and plot PSTH and time course of choice probability
%  using pds and plx files

% subject and experiment session #
subject      = 'nancy';   % monkey name
eList        = 1:37;        % session number. See ephysDataFactory for list.
saveSession  = true;
forceProcess = true;

for expCount = 1:length(eList)
    clearvars -except subject eList saveSession forceProcess expCount
    
    %% check if file already exists for this data set
    %savePath      =
    %'/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/'; %OLD STUF
    % savePath      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/twag/'; % GOODPATH
    savePath      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/twag/allTrialPTA'; % GOODPATH
    
    processedFile = [savePath subject 'Session' num2str(eList(expCount))];
    
    if exist([processedFile '.mat'], 'file') && forceProcess == false;
        
        disp('This file already exists. Loading...')
        load([processedFile '.mat'])
        disp('Done')
    else
        disp('processing... PROCESSING...')
        % call ephysDataFactory
        %[pdsfile, plxfile, PDS, ppNum, PL2PTBfit, PL2PTB, PTB2PL, maxreconstructionerror, flagPlexonTrialStartTime] = ephysDataFactory(subject, eList(expCount));
        [pdsfile, plxfile, PDS, ppNum, synClock] = twagTrimmedDataFactory(subject, eList(expCount));
        
        % choose channels you want to load and analyze. Also load and rename some
        % behavioral variables for convenience
        channelNumber = [1:24] ;
        nCh = length(channelNumber);
        %     allTrSpikes = cell(nCh, 1);
        %     allTrChoice = cell(nCh, 1);
        spikeTimes  = cell(nCh, 1);
        
        pp = ppLoader(subject, ppNum);
        pp = pp{1}.pp;
        
        goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
        goodtrial = cell2mat(goodtrial);
        goodtrial(isnan(goodtrial)) = 0;
        
        stimDistNum = cellfun(@(x) x.theseGabs.stimDistNum, PDS.data, 'UniformOutput', false);
        stimDistNum = cell2mat(stimDistNum);
        %stimDistNum = stimDistNum(goodtrial==1);
        
        sumCoh = cellfun(@(x) x.theseGabs.sumCoh, PDS.data, 'UniformOutput', false);
        sumCoh = cell2mat(sumCoh);
        
        % set the centering time based on desired window
        motionStart = cellfun(@(x) x.stimulus.statesStartTime(4), PDS.data, 'UniformOutput', false);
        motionStart = cell2mat(motionStart)';
        
        if size(syncClock.flagPlexonTrialStartTime) == size(motionStart)  % fix for sessions that run longer after turning of recording. Must vet this.
            centerTime = syncClock.flagPlexonTrialStartTime + motionStart;
        else
            centerTime = syncClock.flagPlexonTrialStartTime + motionStart(1:length(syncClock.flagPlexonTrialStartTime), :);
        end
        
        for jj = 1:nCh
            %% GET SPIKES & ALIGN TO MOTION ON
            [n, spikeTimes{jj, 1}] = plx_ts(plxfile, channelNumber(jj), 1);% FILENAME, CHANNEL, UNIT
        end
        
        %% calculate and plot PSTH
        %  for the five stimulus distributions
        [mSpikesByDist, sSpikesByDist, tSpCnt, dprime, dirPref] = twagPSTH(channelNumber, spikeTimes, stimDistNum, centerTime, sumCoh);
        
        %% calculate and plot time course of CP for each unit and then pulse CP for each unit
        %  from all revco trials and then from frozen trials only
        [mCPallRevco, sCPallRevco, mCPfrozen, sCPfrozen, grandCP_allRevco, grandCP_frozen, bc, allTrSpikes, allTrChoice] = twagCP(channelNumber, spikeTimes, pp, centerTime, goodtrial, dirPref);
        
        % create an index of channels with only a certain d'
        %goodChannel = abs(dprime) > 0.5; %%% this is going to be key to dial in (or just find a new way to label 'good channels')
        %%% might have to just make a sessionCP function and manually enter the channels you want.
        goodChannelIx = abs(dprime) > 2;
        goodChannel   = find(goodChannelIx);
        
        %     manualThresholdSwitch = true;
        %     if manualThresholdSwitch == true
        %         goodChannel = 13:20;
        %         goodChannelIx = zeros(24, 1);
        %         goodChannelIx(goodChannel) = 1;
        %     end
        %%
        % params for separating things by pulse --- pretty hardcoded
        pLength = 150;
        pStart = 350;
        pStartBins = [350 500 650 800 950 1100 1250];
        
        % pulse CP
        for kk = 1:length(goodChannel)
            for ii = 1:7
                pulse_mCP{kk}(ii, :) = mCPallRevco{goodChannel(kk)}(pStartBins(ii):(pStartBins(ii)+pLength));
                %pulse_sCP(ii, :) = sCPallRevco{kk}(pStartBins(ii):(pStartBins(ii)+pLength));
            end
            pulseCPm(:,kk)   = mean(pulse_mCP{kk},2);
            pulseCPstd(:,kk) = std(pulse_mCP{kk}');
        end
        
        %% calculate grand mean CP for the whole session and then pulse CP for whole session
        
        session_mCP_allRevco = cell2mat(mCPallRevco(goodChannelIx));
        session_mCP_allRevco = reshape(session_mCP_allRevco, [1700, sum(goodChannelIx)]);
        session_mCP_allRevco = mean(session_mCP_allRevco,2);
        
        session_sCP_allRevco = cell2mat(sCPallRevco(goodChannelIx));
        session_sCP_allRevco = reshape(session_sCP_allRevco, [1700, sum(goodChannelIx)]);
        session_sCP_allRevco = mean(session_sCP_allRevco,2);
        
        figure; hold on
        plot(bc, session_mCP_allRevco,'k',...
            bc, session_mCP_allRevco+session_sCP_allRevco, 'k--',...
            bc, session_mCP_allRevco-session_sCP_allRevco, 'k--');
        xlim([0 1.5]);
        
        % pulse CP
        for ii = 1:7
            sessionPulse_mCP_allRevco(ii, :) = session_mCP_allRevco(pStartBins(ii):(pStartBins(ii)+pLength));
            %sessionPulse_sCP(ii, :) = session_sCP(pStartBins(ii):(pStartBins(ii)+pLength));
        end
        
        sessionPulseCPm_allRevco   = mean(sessionPulse_mCP_allRevco,2);
        sessionPulseCPstd_allRevco = std(sessionPulse_mCP_allRevco');
        sessionPulseCPse_allRevco   = sessionPulseCPstd_allRevco/sqrt(length(sessionPulse_mCP_allRevco));
        
        %% do it again only using frozen only CP
        
        session_mCP_frozenOnly = cell2mat(mCPfrozen(goodChannelIx));
        session_mCP_frozenOnly = reshape(session_mCP_frozenOnly, [1700, sum(goodChannelIx)]);
        session_mCP_frozenOnly = mean(session_mCP_frozenOnly,2);
        
        session_sCP_frozenOnly = cell2mat(sCPfrozen(goodChannelIx));
        session_sCP_frozenOnly = reshape(session_sCP_frozenOnly, [1700, sum(goodChannelIx)]);
        session_sCP_frozenOnly = mean(session_sCP_frozenOnly,2);
        
        figure; hold on
        plot(bc, session_mCP_frozenOnly,'k',...
            bc, session_mCP_frozenOnly+session_sCP_frozenOnly, 'k--',...
            bc, session_mCP_frozenOnly-session_sCP_frozenOnly, 'k--');
        xlim([0 1.5]);
        
        % pulse CP
        for ii = 1:7
            sessionPulse_mCP_frozenOnly(ii, :) = session_mCP_frozenOnly(pStartBins(ii):(pStartBins(ii)+pLength));
            %sessionPulse_sCP(ii, :) = session_sCP(pStartBins(ii):(pStartBins(ii)+pLength));
        end
        
        sessionPulseCPm_frozenOnly = mean(sessionPulse_mCP_frozenOnly,2);
        sessionPulseCPstd_frozenOnly = std(sessionPulse_mCP_frozenOnly');
        sessionPulseCPse_frozenOnly = sessionPulseCPstd_frozenOnly/sqrt(length(sessionPulse_mCP_frozenOnly));
        
        %% Pulse triggered average
        if length(pp.goodtrial) == length(centerTime)
            centerTimeGood = centerTime(pp.goodtrial==1);
            centerTimeGood = centerTimeGood(pp.stimDistNum==3);
        else
            goodtrial   = pp.goodtrial(1:length(centerTime));
            stimDistNum = stimDistNum(1:length(centerTime));
            stimDistNum = stimDistNum(goodtrial==1);
            
            centerTimeGood = centerTime(goodtrial==1);
            centerTimeGood = centerTimeGood(stimDistNum==3);
        end
        
        [PTA, maxPTA] = pulseTriggeredAverage(pp.cohmat, pp.stimDistNum, spikeTimes, centerTimeGood, goodChannel, dirPref);
        maxPTAm = mean(maxPTA);
        
        %% summary plot using all revco
        %plot kernel and pulse cp
        figure
        subplot(2,3,1);
        errorbar(1:7, pp.kernel.b, pp.kernel.s.se, ...
            'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
        title('Psychophysical kernel')
        xlabel('Pulse')
        ylabel('Weight')
        set(gca, 'Xtick', 1:7)
        xlim([.5 7.5]);
        
        subplot(2,3,2);
        errorbar(1:7, sessionPulseCPm_allRevco, sessionPulseCPse_allRevco, ...
            'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
        title('Choice probability by Pulse')
        h = refline(0,0.5);
        set(h, 'Color', [0 0 0], 'LineStyle', '--');
        xlabel('Pulse')
        ylabel('CP')
        set(gca, 'Xtick', 1:7)
        xlim([.5 7.5]);
        
        subplot(2,3,3);
        plot(1:7, maxPTAm,...
            'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
        title('Max PTA value by pulse')
        xlabel('Pulse')
        ylabel('Pta max value')
        set(gca, 'Xtick', 1:7)
        xlim([.5 7.5]);
        
        % plot them against each other
        subplot(2,3,4);
        scatter(pp.kernel.b, sessionPulseCPm_allRevco);
        xlabel('Behavioral weight');
        ylabel('Pulse CP');
        title('PPK v. CP')
        
        % plot them against each other
        subplot(2,3,5);
        scatter(pp.kernel.b, maxPTAm);
        xlabel('Behavioral weight');
        ylabel('PTA max value');
        title('PPK v. PTA')
        
        % plot them against each other
        subplot(2,3,6);
        scatter(sessionPulseCPm_allRevco, maxPTAm);
        xlabel('Pulse CP');
        ylabel('PTA max value');
        title('CP v. PTA')
        
        %% summary plot using frozen only
        %plot kernel and pulse cp
        figure
        subplot(2,3,1);
        errorbar(1:7, pp.kernel.b, pp.kernel.s.se, ...
            'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
        title('Psychophysical kernel')
        xlabel('Pulse')
        ylabel('Weight')
        set(gca, 'Xtick', 1:7)
        xlim([.5 7.5]);
        
        subplot(2,3,2);
        errorbar(1:7, sessionPulseCPm_frozenOnly, sessionPulseCPse_frozenOnly, ...
            'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
        title('Choice probability by Pulse')
        h = refline(0,0.5);
        set(h, 'Color', [0 0 0], 'LineStyle', '--');
        xlabel('Pulse')
        ylabel('CP')
        set(gca, 'Xtick', 1:7)
        xlim([.5 7.5]);
        
        subplot(2,3,3);
        plot(1:7, maxPTAm,...
            'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
        title('Max PTA value by pulse')
        xlabel('Pulse')
        ylabel('Pta max value')
        set(gca, 'Xtick', 1:7)
        xlim([.5 7.5]);
        
        % plot them against each other
        subplot(2,3,4);
        scatter(pp.kernel.b, sessionPulseCPm_frozenOnly);
        xlabel('Behavioral weight');
        ylabel('Pulse CP');
        title('PPK v. CP')
        
        % plot them against each other
        subplot(2,3,5);
        scatter(pp.kernel.b, maxPTAm);
        xlabel('Behavioral weight');
        ylabel('PTA max value');
        title('PPK v. PTA')
        
        % plot them against each other
        subplot(2,3,6);
        scatter(sessionPulseCPm_frozenOnly, maxPTAm);
        xlabel('Pulse CP');
        ylabel('PTA max value');
        title('CP v. PTA')
        
        
        %% save it
        if saveSession == true
            disp('Saving ephys session...');
            save([savePath subject 'Session' num2str(eList(expCount))]);
            disp('Session saved');
        end
        
    end
end
