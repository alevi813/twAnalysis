%% ephys data analysis for twag
%  using pds and plx files
clear all
% subject and session #
subject      = 'nancy';   % monkey name
eList        = 13;         % session number. See ephysDataFactory for list.

figDir = '/Users/Aaron/Dropbox/twagAnalysis4.1/figures';

% set up a bunch of flags for how you would like to analyze/view the data
saveSession  = false; % save?
forceProcess = true; % process this if a file already exists?
psth         = true; % view the psth for all channels?
cpPlot       = true; % plot individual ch CP time course?
ptaPlot      = true; % plot the PTA?
%%
for expCount = 1:length(eList)
    % check if file already exists for this data set
    % using the thresholded session set
    savePath      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/twag/'; %normal dir 
    %savePath      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/twag/allTrialPTA'; % test dir for pta with all trials
    processedFile = [savePath subject 'Session' num2str(eList(expCount))];
    
    if exist([processedFile '.mat'], 'file') && forceProcess == false
        
        disp('This file already exists. Loading...')
        load([processedFile '.mat'])
        disp('Done')
    else
        disp('processing... PROCESSING...')
        % call dataFactory -- using 'trimmed' version (thresholded session set)
        %[pdsfile, plxfile, PDS, ppNum, ~, ~, ~, ~, flagPlexonTrialStartTime] = twagTrimmedDataFactory(subject, eList(expCount));
        [pdsfile, plxfile, PDS, ppNum, syncClock] = twagTrimmedDataFactory(subject, eList(expCount));

        
        % Load all channels. Rename some behavioral variables for
        % convenience. Preallocate spikeTimes.
        if strcmp(subject, 'nancy')        
            channelNumber = 1:24;
        else
            channelNumber = 1:32;
        end
        nCh = length(channelNumber);
        spikeTimes  = cell(nCh, 1);
        
        pp = ppLoader(subject, ppNum);
        pp = pp{1}.pp;
        
        goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
        goodtrial = cell2mat(goodtrial);
        goodtrial(isnan(goodtrial)) = 0;
        
        stimDistNum = cellfun(@(x) x.theseGabs.stimDistNum, PDS.data, 'UniformOutput', false);
        stimDistNum = cell2mat(stimDistNum);        
        
        sumCoh = cellfun(@(x) x.theseGabs.sumCoh, PDS.data, 'UniformOutput', false);
        sumCoh = cell2mat(sumCoh);
        
        % set the centering time based on motion start
        motionStart = cellfun(@(x) x.stimulus.statesStartTime(4), PDS.data, 'UniformOutput', false);
        motionStart = cell2mat(motionStart)';
        
        % fix for sessions that run longer after turning of recording.
        if size(syncClock.plexonTrialStartTime) == size(motionStart)
            centerTime = syncClock.plexonTrialStartTime + motionStart;
        else
            centerTime = syncClock.plexonTrialStartTime + motionStart(1:length(syncClock.plexonTrialStartTime), :);
        end
        
        %% GET SPIKES & ALIGN TO MOTION ON
        for jj = 1:nCh
%            [n, spikeTimes{jj, 1}] = plx_ts(plxfile, channelNumber(jj), 1);% FILENAME, CHANNEL, UNIT
            [n, spikeTimes{jj, 1}] = plx_ts(plxfile, channelNumber(jj), 0);% FILENAME, CHANNEL, UNIT
        end
        
        %% calculate and plot PSTH for five stimulus distributions
        if psth == true
            [mSpikesByDist, sSpikesByDist, tSpCnt, dprime, dirPref] = twagPSTH(1, channelNumber, spikeTimes, stimDistNum, centerTime, sumCoh);
            set(gcf, 'PaperSize', [12 16], 'PaperPosition', [0 0 12 16])
            %saveas(gcf, [figDir '/psth/psth_channel_session' num2str(eList(expCount))], 'epsc');
            close all
        else
            [~, ~, ~, dprime, dirPref] = twagPSTH(0, channelNumber, spikeTimes, stimDistNum, centerTime, sumCoh);
        end
        
        % create an index of channels with only a certain d'
        goodChannelIx = abs(dprime) > 2;
        goodChannel   = find(goodChannelIx);
        
        %% calculate and plot time course of CP for each unit,
        % pulse CP for each unit
        % using all revco trials or frozen trials only        
        [cp, bc, allTrSpikes, allTrChoice] = twagChannelCP(cpPlot, channelNumber, spikeTimes, pp, centerTime, goodtrial, dirPref);
        if cpPlot            
            set(figure(1), 'PaperSize', [12 16], 'PaperPosition', [0 0 12 16])
            saveas(figure(1), [figDir '/cp/cp_allrevco_channel_session' num2str(eList(expCount))], 'epsc');
            
            set(figure(2), 'PaperSize', [12 16], 'PaperPosition', [0 0 12 16])
            %saveas(figure(2), [figDir '/cp/cp_frozen_chanell_session' num2str(eList(expCount))], 'epsc');
            
            close all
        end
        
        % params for separating things by pulse --- hardcoded...
        pLength = 150;
        pStart = 350;
        pStartBins = [350 500 650 800 950 1100 1250];
        
        % pulse CP
        for kk = 1:length(goodChannel)
            for ii = 1:7
                pulse_mCP{kk}(ii, :) = cp.channel.allRevco.m{goodChannel(kk)}(pStartBins(ii):(pStartBins(ii)+pLength));
                %pulse_sCP(ii, :) = sCPallRevco{kk}(pStartBins(ii):(pStartBins(ii)+pLength));
            end
            cp.channel.pulse.m(:,kk)   = mean(pulse_mCP{kk},2);
            cp.channel.pulse.std(:,kk) = std(pulse_mCP{kk}');
        end
        
        %% session CP
        cp = twagSessionCP(cpPlot, cp, goodChannelIx, bc);
        if cpPlot
            set(figure(1), 'PaperSize', [6 6], 'PaperPosition', [0 0 6 6])
            saveas(figure(1), [figDir '/cp/cp_allrevco_avg_session' num2str(eList(expCount))], 'epsc');
            
            set(figure(2), 'PaperSize', [6 6], 'PaperPosition', [0 0 6 6])
            %saveas(figure(2), [figDir '/cp/cp_frozen_avg_session' num2str(eList(expCount))], 'epsc');
            
            close all
        end
        
        % pulse CP all revCo
        for ii = 1:7
            sessionPulse_mCP_allRevco(ii, :) = cp.session.allRevco.m(pStartBins(ii):(pStartBins(ii)+pLength));
            %sessionPulse_sCP(ii, :) = session_sCP(pStartBins(ii):(pStartBins(ii)+pLength));
        end
        
        cp.session.allRevco.pulse.m    = mean(sessionPulse_mCP_allRevco,2);
        cp.session.allRevco.pulse.std  = std(sessionPulse_mCP_allRevco');
        cp.session.allRevco.pulse.se   = cp.session.allRevco.pulse.std/sqrt(length(sessionPulse_mCP_allRevco));
        
        % pulse CP frozen
        for ii = 1:7
            sessionPulse_mCP_frozen(ii, :) = cp.session.frozen.m(pStartBins(ii):(pStartBins(ii)+pLength));
            %sessionPulse_sCP(ii, :) = session_sCP(pStartBins(ii):(pStartBins(ii)+pLength));
        end
        
        cp.session.frozen.pulse.m    = mean(sessionPulse_mCP_frozen,2);
        cp.session.frozen.pulse.std  = std(sessionPulse_mCP_frozen');
        cp.session.frozen.pulse.se   = cp.session.frozen.pulse.std/sqrt(length(sessionPulse_mCP_frozen));
        %% Pulse triggered average
        if length(pp.goodtrial) == length(centerTime)       % fix for sessions that run longer after turning of recording.
            centerTimeGood = centerTime(pp.goodtrial==1);
            %centerTimeGood = centerTimeGood(pp.stimDistNum==3);
        else
            goodtrial   = pp.goodtrial(1:length(centerTime));
            stimDistNum = stimDistNum(1:length(centerTime));
            stimDistNum = stimDistNum(goodtrial==1);
            
            centerTimeGood = centerTime(goodtrial==1);
            %centerTimeGood = centerTimeGood(stimDistNum==3);
        end
        
        % be careful of all trial vs rc only confusion w/respect to
        % centerTime and cohmat
        [PTA, maxPTA] = pulseTriggeredAverage(ptaPlot, pp.cohmat, pp.stimDistNum, spikeTimes, centerTimeGood, goodChannel, dirPref);
        % [PTA, maxPTA] = pulseTriggeredAverage(ptaPlot, pp.cohmat, pp.stimDistNum, spikeTimes, centerTimeGood(stimDistNum==3), goodChannel, dirPref);
        maxPTAm = mean(maxPTA);
        
        if ptaPlot
            set(figure(1), 'PaperSize', [16 12], 'PaperPosition', [0 0 16 12])
            %saveas(figure(1), [figDir '/ptaAllTrial/pta_channel_session' num2str(eList(expCount))], 'epsc');

            close all
        end
        %% summary plots
        %twagSummaryPlot(pp, cp, maxPTAm);
        close all
        
        %% save it
        if saveSession == true
            disp('Saving ephys session...');
            save([savePath subject 'Session' num2str(eList(expCount))]);
            disp('Session saved');
        end
    end
end
