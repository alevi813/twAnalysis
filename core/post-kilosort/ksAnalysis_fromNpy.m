
savefigs = 1;
postmanualsort = 1;
ptaPlot = 0;

fullDir = ('/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/nancy');
sessions = dir(fullDir);
sessions=sessions(~ismember({sessions.name},{'.','..','.DS_Store','allSessions.mat'})); % CLEAN DAT SHIT UP, also keep it from crashing if they already ran it and have a sessions.mat

sessions = sessions(9); % testing

%%
for nSession = 1:length(sessions)
    dateDir = ([fullDir filesep sessions(nSession).name]);
    expts = dir(dateDir);
    expts = expts(~ismember({expts.name},{'.','..','.DS_Store'})); % CLEAN UP
    
    if length(expts) > 1
        expts = expts(1);
    end
    
    for nExpt = 1:length(expts)
        %% load  stuff
        %  PDS, clock info, & rez
        exDir = ([dateDir filesep expts(nExpt).name]);
        
        load([exDir filesep 'PDS.mat']);
        load([exDir filesep 'syncClock.mat']);
        load([exDir filesep 'KiloSort' filesep 'rez.mat']);
        
        % annoying workaround for loading pp...
        f = dir(exDir);
        f = f(~ismember({f.name},{'.','..','.DS_Store'})); % CLEAN UP
        
        for ii = 1:length(f)
            f(ii).ispp = strcmp(f(ii).name(end-5:end), 'pp.mat');
        end
        ispp = arrayfun(@(x) x.ispp, f, 'uniformoutput', 0);
        ispp = cell2mat(ispp);
        f = f(ispp);
        load([exDir filesep f.name]);
        
        % set paths to spike times and clusters, load em.
        npy.spikeTimes    = ([exDir filesep 'KiloSort' filesep 'spike_times.npy']);
        npy.spikeClusters = ([exDir filesep 'KiloSort' filesep 'spike_clusters.npy']);
        
        %% get spike times and clusters from kilo.rez
        sampleRate = rez.ops.maxFR;
        spikeTimes = readNPY(npy.spikeTimes);
        spikeTimes = double(spikeTimes);
        spikeTimes = spikeTimes./sampleRate;
        
        spikeTimesPtb = nsync.PL2PTB(spikeTimes);
        
        clusters   =  readNPY(npy.spikeClusters);
        
        % read in cluster info from the csv file with group labels, if it exists.
        if exist([exDir filesep 'KiloSort' filesep 'cluster_groups.csv'], 'file')
            [clusterID, clusterGroup] = readClusterCSV(exDir);
            clusterGood = clusterID(strcmp(clusterGroup, 'good'));
            clusterMua  = clusterID(strcmp(clusterGroup, 'mua'));
            
            clusterID = sort([clusterGood; clusterMua]);
        else
            clusterID  = unique(clusters);
        end
        nClust = length(clusterID);
        
        
        % get necessary vars from PDS
        goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
        goodtrial = cell2mat(goodtrial);
        goodtrial(isnan(goodtrial)) = 0;
        
        stimDistNum = cellfun(@(x) x.theseGabs.stimDistNum, PDS.data, 'UniformOutput', false);
        stimDistNum = cell2mat(stimDistNum);
        
        sumCoh = cellfun(@(x) x.theseGabs.sumCoh, PDS.data, 'UniformOutput', false);
        sumCoh = cell2mat(sumCoh);
        
        % set the centering time based on motion start
        motionStart = cellfun(@(x) x.stimulus.statesStartTime(4),  PDS.data, 'UniformOutput', false);
        motionStart = cell2mat(motionStart)';
        
        if size(nsync.computerTrialStartTime) == size(motionStart)
            centerTime = nsync.computerTrialStartTime + motionStart;
        else
            centerTime = nsync.computerTrialStartTime + motionStart(1:length(nsync.computerTrialStartTime), :);
        end
        
        
        %% psth by motion distribution (let's make this its own function!!!)
        figure;
        nrows = nClust/4;
        nrows = floor(nrows) + rem(nClust, 4);
        ixix = [1 2 3 4 5];
        clr  = [1 0 0; 1 .6 .6; 0 0 0; .6 .6 1; 0 0 1];
        for ii = 1:nClust
            theseSpikes = double(spikeTimesPtb(clusters == clusterID(ii)));
            
            for kk = 1:length(ixix)
                distIx   = stimDistNum(:)==ixix(kk);
                if length(distIx) ~= length(centerTime)     % fix for sessions that run longer after turning of recording
                    distIx = distIx(1:length(centerTime));
                end
                
                if sum(theseSpikes) > 0 % whhhhhyyyyyyyyY??
                    [m,s,bc,~,~] = pdsa.eventPsth(theseSpikes, centerTime(distIx), [-.2 1.5],  .001, ones(100,1)/100);
                end
                
                hold on
                subplot(nrows, 4, ii)
                plot(bc, m, 'Color', clr(kk,:), 'LineWidth', 1.75)
                xlim([-.2 1.5])
                supertitle(num2str(dateDir(end-7:end)), 12);
                
            end
            title(['ID ' num2str(clusterID(ii))]);
            
            % CALCULATE d' (let's make this it's own function too!!!)
            ixLeft  = sumCoh < 0;
            ixRight = sumCoh > 0;
            if length(ixLeft) ~= length(centerTime)     % fix for sessions that run longer after turning of recording
                ixLeft  = ixLeft(1:length(centerTime));
                ixRight = ixRight(1:length(centerTime));
            end
            
            if sum(theseSpikes) > 0 % whhhhhyyyyyyyyY??
                [mLeftTrials,sLeftTrials,~, ~, ~]=pdsa.eventPsth(theseSpikes, centerTime(ixLeft), [-.2 1.5], .001, ones(100,1)/100);
                [mRightTrials,sRightTrials, ~, ~, ~]=pdsa.eventPsth(theseSpikes, centerTime(ixRight), [-.2 1.5], .001, ones(100,1)/100);
            end
            
            mL = mean(mLeftTrials); %strong left
            sL = mean(sLeftTrials);
            mR = mean(mRightTrials); %strong right
            sR = mean(sRightTrials);
            
            seDP = sqrt((sL^2 + sR^2)/2);
            
            dprime{nSession}(ii) = (mR - mL)/seDP;
            dprime{nSession}(ii) = round(dprime{nSession}(ii), 2);
            % ^^^ Let's set it up so that:
            % positive dprime means right directional preference
            % negative dprime means left directional preference
            if dprime{nSession}(ii) > 0
                dirPref{nSession}(ii) = 1; %right
            else
                dirPref{nSession}(ii) = 0; %left
            end
            
        end % cluster loop for plotting by dist
        
        %% calculate and plot cp time course
        if length(goodtrial)==length(centerTime)     % fix for sessions that run longer after turning of recording. Must vet this.
            centerTimeGT = centerTime(goodtrial==1);
        else
            goodtrial = goodtrial(1:length(centerTime));
            centerTimeGT = centerTime(goodtrial==1);
        end
        
        figure;
        % assign choice indices according to directional preference
        for ii = 1:nClust
            if dirPref{nSession}(ii) == 1
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
            
            % for ii = 1:nClust
            theseSpikes = double(spikeTimesPtb(clusters == clusterID(ii)));
            
            % align
            [mT1,sT1,bc, ~, trSpkCntT1]=pdsa.eventPsth(theseSpikes, centerTimeGT(T1ix), [-.2 1.5], .001, ones(100,1)/100);
            [mT2,sT2, ~, ~, trSpkCntT2]=pdsa.eventPsth(theseSpikes, centerTimeGT(T2ix), [-.2 1.5], .001, ones(100,1)/100);
            
            allTrSpikes = [trSpkCntT1; trSpkCntT2];
            allTrChoice = [ones(size(trSpkCntT1, 1),1); zeros(size(trSpkCntT2, 1),1)];
            
            % calculate cp
            [mCPallRevco(ii,:), sCPallRevco(ii,:)] = choiceProbabilityCalculate(allTrSpikes, allTrChoice);
            
            % plot
            subplot(nrows, 4, ii)
            hold on
            plot(bc, mCPallRevco(ii,:),'LineWidth', 1.75, 'Color', [0 0 0])
            plot(bc, (mCPallRevco(ii,:)+sCPallRevco(ii,:)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
            plot(bc, (mCPallRevco(ii,:)-sCPallRevco(ii,:)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
            
            title(['ID ' num2str(clusterID(ii))]);
            h = refline(0, 0.5);
            set(h, 'color', 'r')
            xlim([0 1.5]);
            supertitle(num2str(dateDir(end-7:end)), 12);
            
            %%
            bad = num2str([9   6]);
            
            if abs(dprime{nSession}(ii)) > 1 && ~strcmp(bad, [nSession, ii])
                
                [PTA, maxPTA] = pta_postKilo(ptaPlot, pp.cohmat, pp.stimDistNum, theseSpikes, centerTimeGT, dirPref{nSession}(ii));
                %maxPTAm = mean(maxPTA,2);
            else
                PTA = nan;
                maxPTA = nan;
            end
            allPTA{nSession, ii} = PTA;
            allMaxPTA{nSession, ii} = maxPTA;
            
            
        end % cluster loop for calc and plotting
        
    end % experiment (nExpts, within dateDir) loop
    
    goodClustIx{nSession} = abs(dprime{nSession}) > 1;
    if sum(goodClustIx{nSession}) > 1
        mSessionCP_allRevco(nSession, :)    = mean(mCPallRevco(goodClustIx{nSession}, :));
        %        sSessionCP_allRevco(nSession, :)     = std(mCPallRevco(goodClustIx{nSession}, :));
        sSessionCP_allRevco(nSession, :)     = mean(sCPallRevco(goodClustIx{nSession}, :));
    else
        mSessionCP_allRevco(nSession, :)    = mCPallRevco(goodClustIx{nSession}, :);
        sSessionCP_allRevco(nSession, :)    = sCPallRevco(goodClustIx{nSession}, :);
    end
    
    figure
    hold on
    plot(bc, mSessionCP_allRevco(nSession, :),'LineWidth', 1.75, 'Color', [0 0 0])
    plot(bc, (mSessionCP_allRevco(nSession, :)+  sSessionCP_allRevco(nSession, :)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
    plot(bc, (mSessionCP_allRevco(nSession, :)-  sSessionCP_allRevco(nSession, :)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
    title(['ID ' num2str(clusterID(ii))]);
    h = refline(0, 0.5);
    set(h, 'color', 'r')
    
    
    
    % params for separating things by pulse --- hardcoded...
    pLength = 150;
    pStart = 350;
    pStartBins = [350 500 650 800 950 1100 1250];
    
    % pulse CP
    for ii = 1:7
        pulse_mCP{nSession}(ii, :) = mSessionCP_allRevco(pStartBins(ii):(pStartBins(ii)+pLength));
        %pulse_sCP(ii, :) = conditionMean(pStartBins(ii):(pStartBins(ii)+pLength));
    end
    mPulseCP(:, nSession)  = mean(pulse_mCP{nSession},2);
    sPulseCP(:, nSession)  = std(pulse_mCP{nSession}');
    
end % session (nSessions, dateDir) loop

%%
conditionMean = mean(mSessionCP_allRevco);
conditionStd = mean(sSessionCP_allRevco);

figure; hold on
plot(bc, conditionMean, 'LineWidth', 2, 'Color', [0 0 0 ])
plot(bc, conditionMean + conditionStd, 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
plot(bc, conditionMean - conditionStd, 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);

xlim([0 1.2])

% params for separating things by pulse --- hardcoded...
pLength = 150;
pStart = 350;
pStartBins = [350 500 650 800 950 1100 1250];

% pulse CP
for ii = 1:7
    condition_pulse_mCP(ii, :) = conditionMean(pStartBins(ii):(pStartBins(ii)+pLength));
    %pulse_sCP(ii, :) = conditionMean(pStartBins(ii):(pStartBins(ii)+pLength));
end
mPulseCP  = mean(condition_pulse_mCP,2);
sPulseCP  = std(condition_pulse_mCP');



%% save some figs
%%% hardcoded %%%
% if postmanualsort
%     baseDir = '/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/figures/post manual sort/';
% else
%     baseDir = '/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/figures/';
% end
%
% for ii = 1:(nSession*2)
%     figure(ii)
%
%     if mod(ii,2) == 0
%         saveas(gcf, [baseDir 'cp' num2str(ii)], 'epsc');
%     else
%         saveas(gcf, [baseDir 'psth' num2str(ii)], 'epsc');
%     end
% end

