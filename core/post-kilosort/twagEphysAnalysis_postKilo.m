savedat = false;
savefigs = false;

[~, host] = system('hostname');

if strcmp(host(1), 'd') %dhcp-129-116-178-237.cps.utexas.edu
    %desktop
    figDir = '/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/figures/leo/';
else
    %laptop
    figDir = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/figures/leo/';
end


%fullDir = ('/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/nancy');
%fullDir = ('/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/leo');

subject = 'leo';

if strcmp(subject, 'leo')
    fullDir = ('/Volumes/HukLab/Macaque/Projects/twag/leo');
else
    fullDir = ('/Volumes/HukLab/Macaque/Projects/twag/nancy');
end

sessions = dir(fullDir);
sessions=sessions(~ismember({sessions.name},{'.','..','.DS_Store','allSessions.mat'})); % CLEAN DAT SHIT UP, also keep it from crashing if they already ran it and have a sessions.mat

%sessions = sessions(1:37);
sessionList = 5; % be careful about sessions in your directory with no twagPDS file.
sessions = sessions(sessionList);

%screenLatency = .0524634; %52 ms screen latency for old LG/rig1 set up.
screenLatency = 0;

%%

for nSession = 1:length(sessions)
    %for nSession = 2
    dateDir = ([fullDir filesep sessions(nSession).name]);
    expts = dir(dateDir);
    expts = expts(~ismember({expts.name},{'.','..','.DS_Store'})); % CLEAN UP
    
    %     % there are one or two sessions with more than one expt. Currently hardcoding to only look at first right now.
    %     if length(expts) > 1
    %         expts = expts(1);
    %     end
    
    for ii = 1:length(expts)
        %isFinal(ii) = strcmp(expts(ii).name(end-4:end), 'final');
        isFinal(ii) = strcmp(expts(ii).name(end-3:end), 'test');
    end
    
    if size(expts, 1) > 1
        expts = expts(isFinal);
    end
    
    for nExpt = 1:length(expts)
        %% load  stuff
        %  PDS, clock info, & rez
        exDir = ([dateDir filesep expts(nExpt).name]);
        
        %load([exDir filesep 'PDS.mat']);
        load([exDir filesep 'syncClock.mat']);
        load([exDir filesep 'KiloSort' filesep 'rez.mat']);
        %        load([exDir filesep 'KiloSort2' filesep 'rez2.mat']);
        
        % annoying workaround for loading pp...
        f = dir(exDir);
        f = f(~ismember({f.name},{'.','..','.DS_Store'})); % CLEAN UP
        
        for ii = 1:length(f)
            f(ii).isPDS = strcmp(f(ii).name(end-3:end), '.PDS');
        end
        isPDS = arrayfun(@(x) x.isPDS, f, 'uniformoutput', 0);
        isPDS = cell2mat(isPDS);
        f = f(isPDS);

        % if you have multiple pds files, use the largest one
        for ii = 1:length(f)
            byteSize(ii) = f(ii).bytes;
        end        
        [~, sizeIx] = max(byteSize);
        f = f(sizeIx);
        
        %load([f.folder filesep f.name], '-mat'); % load PDS
        load([exDir filesep f.name], '-mat'); % load PDS
        
        % process PDS
        pp = analyzeBehaviorSession(exDir, f.name, 1, false);
        
        % set paths to spike times and clusters, load em.
        npy.spikeTimes    = ([exDir filesep 'KiloSort' filesep 'spike_times.npy']);
        npy.spikeClusters = ([exDir filesep 'KiloSort' filesep 'spike_clusters.npy']);
        
        %% get spike times and clusters from kilo.rez
        if ~exist('syncClock', 'var')
            syncClock = nsync;
        end
        
        %sampleRate = rez.ops.maxFR;
        sampleRate = rez.ops.fs; % uhhhh 20k for nancy?
        spikeTimes = readNPY(npy.spikeTimes);
        spikeTimes = double(spikeTimes);
        spikeTimes = spikeTimes./sampleRate;
        
        spikeTimesPtb = syncClock.PL2PTB(spikeTimes);
        
        clusters   =  readNPY(npy.spikeClusters);
        
        %read in cluster info from the csv file with group labels, if it exists.
        %if exist([exDir filesep 'KiloSort' filesep 'cluster_groups.csv'], 'file')
        %    [clusterID, clusterGroup] = readClusterCSV(exDir);
        %    clusterGood = clusterID(strcmp(clusterGroup, 'good'));
        %    clusterMua  = clusterID(strcmp(clusterGroup, 'mua'));
        %    
        %    clusterID = sort([clusterGood; clusterMua]);
        %else
            clusterID  = unique(clusters);
        %end
        nClust = length(clusterID);
        
        
        % get necessary vars from PDS
        %         goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
        %         goodtrial = cell2mat(goodtrial);
        %         goodtrial(isnan(goodtrial)) = 0;
        goodtrial = pp.goodtrial;
        
        stimDistNum = cellfun(@(x) x.theseGabs.stimDistNum, PDS.data, 'UniformOutput', false);
        stimDistNum = cell2mat(stimDistNum);
        
        sumCoh = cellfun(@(x) x.theseGabs.sumCoh, PDS.data, 'UniformOutput', false);
        sumCoh = cell2mat(sumCoh);
        
        % set the centering time based on motion start
        motionStart = cellfun(@(x) x.stimulus.statesStartTime(4),  PDS.data, 'UniformOutput', false);
        motionStart = cell2mat(motionStart)';
        
        if size(syncClock.computerTrialStartTime) == size(motionStart)
            centerTime = syncClock.computerTrialStartTime + motionStart;
        else
            centerTime = syncClock.computerTrialStartTime + motionStart(1:length(syncClock.computerTrialStartTime), :);
        end
        
        centerTime = centerTime + screenLatency; % add screen latency if this data was collected on rig 1 with old lg tv.
        
        %% psth by motion distribution (let's make this its own function!!!)
        figure;
        nrows = nClust/4;
        nrows = floor(nrows) + rem(nClust, 4);
        ixix = [1 2 3 4 5];
        clr  = [1 0 0; 1 .6 .6; 0 0 0; .6 .6 1; 0 0 1];
        for ii = 1:nClust
            %figure
            theseSpikes = double(spikeTimesPtb(clusters == clusterID(ii)));
            maxrate(ii) = 0 ; % preset max rate to be filled at end of dist plot loop
            
            for kk = 1:length(ixix)
                distIx   = stimDistNum(:)==ixix(kk) & goodtrial(:) == 1;
                if length(distIx) ~= length(centerTime)     % fix for sessions that run longer after turning off recording
                    distIx = distIx(1:length(centerTime));
                end
                
                if sum(theseSpikes) > 0 % whhhhhyyyyyyyyY??
                    [m,s,bc,v,~] = pdsa.eventPsth(theseSpikes, centerTime(distIx), [-.5 1.5],  .001, ones(100,1)/100);
                    %[m,s,bc,~,~] = pdsa.eventPsth(theseSpikes, centerTime(distIx), [-.2 1.5],  .001, ones(10,1)/10);
                end
                
                hold on
                subplot(nrows, 4, ii);
                plot(bc, m, 'Color', clr(kk,:), 'LineWidth', 1.75);
                xlim([-.25 1.5])
                %supertitle(num2str(dateDir(end-7:end)), 12);
                
                if maxrate(ii) < max(m)
                    maxrate(ii) = max(m);
                end
            end
            %title(['ID ' num2str(clusterID(ii))]);
            box off
            
            %             % CALCULATE d' (let's make this it's own function too!!!)
            %             ixLeft  = sumCoh < 0;
            %             ixRight = sumCoh > 0;
            %             if length(ixLeft) ~= length(centerTime)     % fix for sessions that run longer after turning of recording
            %                 ixLeft  = ixLeft(1:length(centerTime));
            %                 ixRight = ixRight(1:length(centerTime));
            %             end
            %
            %             if sum(theseSpikes) > 0 % whhhhhyyyyyyyyY??
            %                 [mLeftTrials,sLeftTrials,~, ~, ~]=pdsa.eventPsth(theseSpikes, centerTime(ixLeft), [-.5 1.5], .001, ones(100,1)/100);
            %                 [mRightTrials,sRightTrials, ~, ~, ~]=pdsa.eventPsth(theseSpikes, centerTime(ixRight), [-.5 1.5], .001, ones(100,1)/100);
            %             end
            %
            %             if strcmp(PDS.initialParametersMerged.gabs.condition, 'baseline')
            %                 mL = mean(mLeftTrials); %strong left
            %                 sL = mean(sLeftTrials);
            %                 mR = mean(mRightTrials); %strong right
            %                 sR = mean(sRightTrials);
            %             elseif strcmp(PDS.initialParametersMerged.gabs.condition, 'late')
            %                 mL = mean(mLeftTrials(1, 1100:1600)); %strong left
            %                 sL = mean(sLeftTrials(1, 1100:1600));
            %                 mR = mean(mRightTrials(1, 1100:1600)); %strong right
            %                 sR = mean(sRightTrials(1, 1100:1600));
            %             else strcmp(PDS.initialParametersMerged.gabs.condition, 'early')
            %                 mL = mean(mLeftTrials(1, 500:1000)); %strong left
            %                 sL = mean(sLeftTrials(1, 500:1000));
            %                 mR = mean(mRightTrials(1, 500:1000)); %strong right
            %                 sR = mean(sRightTrials(1, 500:1000));
            %             end
            %             seDP = sqrt((sL^2 + sR^2)/2);
            %
            %             dprime(ii) = (mR - mL)/seDP;
            %             dprime(ii) = round(dprime(ii), 2);
            %             % ^^^ Let's set it up so that:
            %             % positive dprime means right directional preference
            %             % negative dprime means left directional preference
            %             if dprime(ii) > 0
            %                 dirPref(ii) = 1; %right
            %             else
            %                 dirPref(ii) = 0; %left
            %             end
            
            [dprime(ii), dirPref(ii)] = twag_dprime(theseSpikes, sumCoh, centerTime, PDS.initialParametersMerged.gabs.condition);
            %[dprime(ii), dirPref(ii)] = twag_dprimeTEST(theseSpikes, sumCoh, centerTime, PDS.initialParametersMerged.gabs.condition);
            
        end % cluster loop for plotting by dist
        
        %% calculate and plot cp time course
        if length(goodtrial)==length(centerTime)     % fix for sessions that run longer after turning of recording. Must vet this.
            centerTimeGT = centerTime(goodtrial==1);
        else
            goodtrial = goodtrial(1:length(centerTime));
            centerTimeGT = centerTime(goodtrial==1);
        end
        
        figure(2); figure(3);
        % assign choice indices according to directional preference
        for ii = 1:nClust
            % get spikes for specific cluster/unit
            theseSpikes = double(spikeTimesPtb(clusters == clusterID(ii)));
            
            % calculate cp!
            [allrevco, frozenonly] = unitCP_postKilo(centerTimeGT, pp, theseSpikes, dirPref(ii));
            % reassign outputs from unitCP_postKilo. annoying but eases use
            cp.allrevco.m(ii,:)   = allrevco.m;
            cp.allrevco.s(ii,:)   = allrevco.s;
            cp.frozenonly.m(ii,:) = frozenonly.m;
            cp.frozenonly.s(ii,:) = frozenonly.s;
            
            % plot allrevco
            figure(2)
            subplot(nrows, 4, ii)
            hold on
            plot(bc,cp.allrevco.m(ii,:),'LineWidth', 1.75, 'Color', [0 0 0])
            plot(bc, (cp.allrevco.m(ii,:)+cp.allrevco.s(ii,:)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
            plot(bc, (cp.allrevco.m(ii,:)-cp.allrevco.s(ii,:)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
            
            title(['ID ' num2str(clusterID(ii))]);
            h = refline(0, 0.5);
            set(h, 'color', 'r')
            xlim([-.5 1.5]);
            supertitle([num2str(dateDir(end-7:end)) ' allrevco'], 12);
            
            % plot frozenonly
            figure(3)
            subplot(nrows, 4, ii)
            hold on
            plot(bc,cp.frozenonly.m(ii,:),'LineWidth', 1.75, 'Color', [0 0 0])
            plot(bc, (cp.frozenonly.m(ii,:)+cp.frozenonly.s(ii,:)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
            plot(bc, (cp.frozenonly.m(ii,:)-cp.frozenonly.s(ii,:)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
            
            title(['ID ' num2str(clusterID(ii))]);
            h = refline(0, 0.5);
            set(h, 'color', 'r')
            xlim([-.5 1.5]);
            supertitle([num2str(dateDir(end-7:end)) ' frozenonly'], 12);
        end % cluster loop for calc and plotting
        
    end % experiment (nExpts, within dateDir) loop
    
    %% threshold by d' and firing rate -- arb. hardcoded number for now.
    dpThresh = 1; maxrateThresh = 10;
    goodClustIx = abs(dprime) > dpThresh & maxrate > maxrateThresh;
    
    figure
    imagesc(goodClustIx)
    colormap(gray)
    %    title(['dprime > ' num2str(dpThresh) ' & maxrate > ' num2str(maxrateThresh)])
    title('manual thresholding')
    
    if sum(goodClustIx) > 1
        cp.allrevco.session.m    = mean(cp.allrevco.m(goodClustIx, :));
        cp.allrevco.session.s    = mean(cp.allrevco.s(goodClustIx, :));
        
        cp.frozenonly.session.m  = mean(cp.frozenonly.m(goodClustIx, :));
        cp.frozenonly.session.s    = mean(cp.frozenonly.s(goodClustIx, :));
    else
        cp.allrevco.session.m   = cp.allrevco.m(goodClustIx, :);
        cp.allrevco.session.s   = cp.allrevco.s(goodClustIx, :);
        
        cp.frozenonly.session.m   = cp.frozenonly.m(goodClustIx, :);
        cp.frozenonly.session.s   = cp.frozenonly.s(goodClustIx, :);
    end
    
    
    %% calculate full session cp (across units). group and calculate mean CP per pulse
    figure
    hold on
    
    if ~isempty(cp.allrevco.session.m)
        
        subplot(2,1,1); hold on
        plot(bc, cp.allrevco.session.m,'LineWidth', 1.75, 'Color', [0 0 0])
        plot(bc, (cp.allrevco.session.m + cp.allrevco.session.s), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
        plot(bc, (cp.allrevco.session.m - cp.allrevco.session.s), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
        h = refline(0, 0.5);
        set(h, 'color', 'r')
        title('allrevco')
        
        % params for separating things by pulse --- hardcoded...
        pLength = 150;
        %         pStart = 300;
        %         pStartBins = [300 450 600 750 900 1050 1200];
        pStart = 500;
        pStartBins = [500 650 800 950 1100 1250 1400];
        
        % pulse CP
        for ii = 1:7
            cp.allrevco.session.pulse.grouped(ii, :) = cp.allrevco.session.m(pStartBins(ii):(pStartBins(ii)+pLength));
        end
        
        cp.allrevco.session.pulse.m  = mean(cp.allrevco.session.pulse.grouped,2);
        cp.allrevco.session.pulse.s  = std(cp.allrevco.session.pulse.grouped');
    end
    
    if ~isempty(cp.frozenonly.session.m)
        
        subplot(2,1,2); hold on
        plot(bc, cp.frozenonly.session.m,'LineWidth', 1.75, 'Color', [0 0 0])
        plot(bc, (cp.frozenonly.session.m + cp.frozenonly.session.s), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
        plot(bc, (cp.frozenonly.session.m - cp.frozenonly.session.s), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
        h = refline(0, 0.5);
        set(h, 'color', 'r')
        title('frozen only')
        
        % params for separating things by pulse --- hardcoded...
        pLength = 150;
        %         pStart = 300;
        %         pStartBins = [300 450 600 750 900 1050 1200];
        pStart = 500;
        pStartBins = [500 650 800 950 1100 1250 1400];
        
        % pulse CP
        for ii = 1:7
            cp.frozenonly.session.pulse.grouped(ii, :) = cp.frozenonly.session.m(pStartBins(ii):(pStartBins(ii)+pLength));
        end
        
        cp.frozenonly.session.pulse.m  = mean(cp.frozenonly.session.pulse.grouped,2);
        cp.frozenonly.session.pulse.s  = std(cp.frozenonly.session.pulse.grouped');
    end
    
    %% save processedData and figures
    if savedat
        
        if strcmp(host(1), 'd') %dhcp-129-116-178-237.cps.utexas.edu
            %desktop
            save(['/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData/leoSession' num2str(sessionList(nSession))]);
        else
            %laptop
            save(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData/leoSession' num2str(sessionList(nSession))]);
        end
        
        disp(['saved session ' num2str(sessionList(nSession))]);
    end
    
    if savefigs
        figure(1)
        saveas(gcf, [figDir 'psth/psth_unit_session' num2str(sessionList(nSession))], 'epsc');
        
        figure(2)
        saveas(gcf, [figDir 'cp/cp_unit_allrevco_session' num2str(sessionList(nSession))], 'epsc');
        
        figure(3)
        saveas(gcf, [figDir 'cp/cp_unit_frozen_session' num2str(sessionList(nSession))], 'epsc');
        
        figure(4)
        saveas(gcf, [figDir 'thresholding/thresholding_session' num2str(sessionList(nSession))], 'epsc');
        
        figure(5)
        saveas(gcf, [figDir 'cp/cp_avg_session' num2str(sessionList(nSession))], 'epsc');
    end
    close all
    clearvars -except nSession fullDir sessions savedat savefigs figDir screenLatency sessionList host
end % session (nSessions, dateDir) loop

%%
% conditionMean = mean(mSessionCP_allRevco);
% conditionStd = mean(sSessionCP_allRevco);
%
% figure; hold on
% plot(bc, conditionMean, 'LineWidth', 2, 'Color', [0 0 0 ])
% plot(bc, conditionMean + conditionStd, 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
% plot(bc, conditionMean - conditionStd, 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
%
% xlim([0 1.2])
%
% % params for separating things by pulse --- hardcoded...
% pLength = 150;
% pStart = 350;
% pStartBins = [350 500 650 800 950 1100 1250];
%
% % pulse CP
% for ii = 1:7
%     condition_pulse_mCP(ii, :) = conditionMean(pStartBins(ii):(pStartBins(ii)+pLength));
%     %pulse_sCP(ii, :) = conditionMean(pStartBins(ii):(pStartBins(ii)+pLength));
% end
% mPulseCP  = mean(condition_pulse_mCP,2);
% sPulseCP  = std(condition_pulse_mCP');

