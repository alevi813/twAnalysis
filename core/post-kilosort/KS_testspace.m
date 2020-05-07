%% get cell arrays of all KiloSort data, clock info, and PDS files
%  either run loadKiloSortSession, or load output its pre-saved output
process  = 0;
fullpsth = 0;


if process
    allSessions = processEphysData('twag', 'nancy', baseDir, 0, 1, 1);
else
    %load
    %     load('/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/nancy/allSessions.mat');
end

fullDir = ('/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/nancy');
sessions = dir(fullDir);
sessions=sessions(~ismember({sessions.name},{'.','..','.DS_Store','allSessions.mat'})); % CLEAN DAT SHIT UP, also keep it from crashing if they already ran it and have a sessions.mat

%% loop over sessions align spikes and plot motion PSTH for each cluster

for jj = 1:length(kilo)
    if jj >= 2
        clearvars -except jj allSessions allPDS clockSync kilo fullDir sessions fullpsth
    end
    
    dateDir = ([fullDir filesep sessions(jj).name]);
    expts = dir(dateDir);
    expts = expts(~ismember({expts.name},{'.','..','.DS_Store'})); % CLEAN UP
    
    for kk = 1:length(expts)
        exDir = ([dateDir filesep expts(kk).name]);
        
        f = dir(exDir);
        f = f(~ismember({f.name},{'.','..','.DS_Store'})); % CLEAN UP
        
        for ii = 1:length(f)
            f(ii).ispp = strcmp(f(ii).name(end-5:end), 'pp.mat');
        end
        ispp = arrayfun(@(x) x.ispp, f, 'uniformoutput', 0);
        ispp = cell2mat(ispp);
        f = f(ispp);
        load([exDir filesep f.name]);
        
        % get spike times and clusters from kilo.rez
        sampleRate = kilo{jj}.rez.ops.maxFR;
        spikeTimes = kilo{jj}.rez.st3(:,1);
        spikeTimes = spikeTimes./sampleRate;
        
        spikeTimesPtb = clockSync{jj}.PL2PTB(spikeTimes);
        
        clusters   = kilo{jj}.rez.st3(:,5);
        clusterID  = unique(clusters);
        nClust = length(clusterID);
        
        % get necessary vars from PDS
        goodtrial = cellfun(@(x) x.pldaps.goodtrial, allPDS{jj}.data, 'UniformOutput', false);
        goodtrial = cell2mat(goodtrial);
        goodtrial(isnan(goodtrial)) = 0;
        
        stimDistNum = cellfun(@(x) x.theseGabs.stimDistNum, allPDS{jj}.data, 'UniformOutput', false);
        stimDistNum = cell2mat(stimDistNum);
        
        sumCoh = cellfun(@(x) x.theseGabs.sumCoh, allPDS{jj}.data, 'UniformOutput', false);
        sumCoh = cell2mat(sumCoh);
        
        % set the centering time based on motion start
        motionStart = cellfun(@(x) x.stimulus.statesStartTime(4),  allPDS{jj}.data, 'UniformOutput', false);
        motionStart = cell2mat(motionStart)';
        
        if size(clockSync{jj}.computerTrialStartTime) == size(motionStart)
            centerTime = clockSync{jj}.computerTrialStartTime + motionStart;
        else
            centerTime = clockSync{jj}.computerTrialStartTime + motionStart(1:length(clockSync{jj}.computerTrialStartTime), :);
        end
        
        %% align to motion on and plot!
        if fullpsth
            figure;
            supertitle(['Day ' num2str(jj) ' Expt ' num2str(kk)],  21) % need to get real date info...
            
            nrows = nClust/4;
            nrows = floor(nrows) + rem(nClust, 4);
            for ii = 1:nClust;
                theseSpikes = double(spikeTimesPtb(clusters == clusterID(ii)));
                
                if sum(theseSpikes) > 0 % whhhhhyyyyyyyyY??
                    [m,s,bc,~,~] = pdsa.eventPsth(theseSpikes, centerTime, [-.2 1.5],  .001, ones(100,1)/100);
                end
                
                subplot(nrows, 4, ii)
                plot(bc, m, 'LineWidth', 2)
                xlim([-.2 1.5])
            end
        end % all trial psth loop
        
        %% psth by motion distribution (let's make this its own function!!!)
        figure;
        supertitle(['Day ' num2str(jj) ' Expt ' num2str(kk)],  21) % need to get real date info...
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
            end
            
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
            
            dprime{jj}(ii) = (mR - mL)/seDP;
            dprime{jj}(ii) = round(dprime{jj}(ii), 2);
            % ^^^ Let's set it up so that:
            % positive dprime means right directional preference
            % negative dprime means left directional preference
            if dprime{jj}(ii) > 0
                dirPref{jj}(ii) = 1; %right
            else
                dirPref{jj}(ii) = 0; %left
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
            if dirPref{jj}(ii) == 1
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
            [mCPallRevco{jj}(ii,:), sCPallRevco{jj}(ii,:)] = choiceProbabilityCalculate(allTrSpikes, allTrChoice);
            
            % plot
            subplot(nrows, 4, ii)
            hold on
            plot(bc, mCPallRevco{jj}(ii,:),'LineWidth', 1.75, 'Color', [0 0 0])
            plot(bc, (mCPallRevco{jj}(ii,:)+sCPallRevco{jj}(ii,:)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
            plot(bc, (mCPallRevco{jj}(ii,:)-sCPallRevco{jj}(ii,:)), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
            
            xlim([0 1.5]);
        end % cluster loop for calc and plotting
    end % expts loop
end % session loop

%% save some figs
%%% hardcoded %%%

baseDir = '/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/figures/';

for ii = 1:72
    figure(ii)
    
    if mod(ii,2) == 0
        saveas(gcf, [baseDir 'cp' num2str(ii)], 'epsc');
    else
        saveas(gcf, [baseDir 'psth' num2str(ii)], 'epsc');
    end 
end
        
%% let's break things down into new functions
% already have these things, but need to make generalizable versions.
% psth
% d' & dirpref 
% cp time course
%

