% clear all

%%
comp    = getComp;
baseDir = '/Volumes/HukLab/Macaque/Projects/twag';

subject = 'leo';
sList   = 28; % 5 -> good flat-stim session 28 -> the good early-stim session

expts = dir([baseDir filesep subject]);
expts = expts(~ismember({expts.name},{'.','..','.DS_Store', 'test', 'toSort'})); % CLEAN UP
expts = expts(sList);

allMaxTheta = [];

for iS = 1:length(sList)
    %fileList = dir([expts(iS).folder filesep expts(iS).name]);
    fileList = dir([baseDir filesep subject filesep expts(iS).name]);
    fileList = fileList(~ismember({fileList.name},{'.','..','.DS_Store'})); % CLEAN UP
    
    for ii = 1:length(fileList)
        isFinal(ii) = strcmp(fileList(ii).name(end-4:end), 'final');
    end
    
    finalDir = fileList(isFinal);
    
    %    fileList = dir([finalDir.folder filesep finalDir.name]);
    fileList = dir([baseDir filesep subject filesep expts(iS).name filesep finalDir.name]);
    fileList = fileList(~ismember({fileList.name},{'.','..','.DS_Store'})); % CLEAN UP
    for ii = 1:length(fileList)
        isPDS(ii) = strcmp(fileList(ii).name(end-3:end), '.PDS');
        isPLX(ii) = strcmp(fileList(ii).name(end-3:end), '.pl2');
    end
    plxfile = fileList(isPLX);
    plxfile = [plxfile.folder filesep plxfile.name];
    
    pdsList = fileList(isPDS);
    %mapFile = pdsList(end); % just take the last one for now.
    mapFile = pdsList(end-2);
    pdsfile = [mapFile.folder filesep mapFile.name];
    
    %%
    % load pdsfile
    load(pdsfile, '-mat');
    %sync clock
    syncClock = syncPlexonClock(PDS, plxfile);
    
    %% take care of ephys
    dataSource = 'kiloSort'; %kiloSort or plx
    
    if strcmp(dataSource, 'kiloSort')
        nCh = 1;
    else
        nCh = 32; % 24 for nancy data... add some options
    end
    
    for iCh = 1:nCh
        switch dataSource
            case 'kiloSort'
                
                if strcmp(comp, 'laptop')
                    neuronDir       = dir('/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo/neurons');
                else
                    neuronDir       = dir('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo/neurons');
                end
                neuronFileNames = arrayfun(@(x) x.name, neuronDir, 'UniformOutput', false);
                neuronDir       = neuronDir(contains(neuronFileNames, expts(iS).name, 'IgnoreCase', true));
                
                nClust = size(neuronDir, 1);
                             
            case 'plx'
                
                channelNumber = iCh;
                nClust        = 1;
                
                % get spike times
                [~, spikeTimes] = plx_ts(plxfile, channelNumber, 0);
                
        end
        
        %%
        nTrials = length(PDS.data);
        for iTrial = 1:nTrials
            
            if isfield(PDS.data{iTrial}.pldaps, 'goodtrial')
                goodtrial(iTrial) = 1;
            else
                goodtrial(iTrial) = nan;
            end
        end
        
        goodtrial(isnan(goodtrial)) = 0;
        
        % get thetas - directions on each trial
        % thetas = cellfun(@(x) x.stimulus.motion1.thetas, PDS.data, 'UniformOutput', false);
        % thetas = cell2mat(thetas);
        % thetas = reshape(thetas, PDS.initialParametersMerged.stimulus.motion1.nThetas, nTrials);
        
        for iTrial = 1:length(PDS.data)
            if isfield(PDS.data{iTrial}.stimulus.motion1, 'thetas')
                thetas(iTrial) = PDS.data{iTrial}.stimulus.motion1.thetas;
            else
                thetas(iTrial) = NaN;
            end
        end
        
        thetaList    = unique(thetas);
        thetaList(isnan(thetaList)) = [];
        % thetaListRad = deg2rad(thetaList);
        thetaListRad = thetaList * pi/180;
        thetaListRad = [thetaListRad(:); thetaListRad(1)]; % wrap it/make it circular
        nThetas      = numel(unique(thetaList));
        
        % get relative motion-on times
        for itrial = 1:length(PDS.data)
            if isfield(PDS.data{itrial}.stimulus.motion1, 'stimStateOnTime') && length(PDS.data{itrial}.stimulus.motion1.stimStateOnTime) > 2
                trialMotionStartTime(itrial) = PDS.data{itrial}.stimulus.motion1.stimStateOnTime(3);
            else
                trialMotionStartTime(itrial) = nan;
            end
        end
        
        % length of time motion is on
        motionLength = PDS.initialParametersMerged.stimulus.motion1.motionStateStartTime(3);
        
        % length of time  a single direction of motion is presented (note: from initial
        % params, not individual trials_
        thetaLength  = motionLength / PDS.initialParametersMerged.stimulus.motion1.nThetas;
        
        % add plexon trial start time to relative motion-on time
        centerTime = syncClock.plexonTrialStartTime + trialMotionStartTime';
        
        % goodtrials only
        thetas     = thetas(:, goodtrial == 1);
        centerTime = centerTime(goodtrial==1);
        
        allThetaTimes(:,1) = centerTime;
        for iTheta = 2:PDS.initialParametersMerged.stimulus.motion1.nThetas
            allThetaTimes(:, iTheta) = centerTime + thetaLength*(iTheta - 1);
        end
        
        %%
        
        win    = [0.05 .25];
        
        % lets start with just the first motion presentation...
        thetas = thetas(1,:);
        
        for iClust = 1:nClust
            
            switch dataSource
                case 'kiloSort'
                    %                 theseSpikes = double(spikeTimes(clusters == clusterID(iClust)));
                    thisNeuron = load([neuronDir(iClust).folder filesep neuronDir(iClust).name]);
                    theseSpikes = thisNeuron.spikeTimes;
                case 'plx'
                    theseSpikes = spikeTimes;
            end
            
            for iTheta = 1:nThetas
                
                % index for one direction at a time
                dirIx = thetas==thetaList(iTheta);
                
                [mSpikes, sSpike ,bc, ~, totSpk] = pdsa.eventPsth(theseSpikes, centerTime(dirIx), win,  .001, ones(100,1)/100);
                %[mSpikes, sSpike ,bc, ~, totSpk] = pdsa.eventPsth(theseSpikes, allThetaTimes(dirIx'), win,  .001, ones(100,1)/100);
                
                mResponseByDir(iClust, iTheta)   = mean(mSpikes);
                sResponseByDir(iClust, iTheta)   = std(mSpikes);
                semResponseByDir(iClust, iTheta) = sResponseByDir(iClust, iTheta)/sqrt(sum(sum(dirIx,2)));
                
                spikeCounts{iClust, iTheta} = nansum(totSpk,2);
                [~, ~, spikeCounts{iClust, iTheta}] = find(spikeCounts{iClust, iTheta});
                
                mSpikeCount(iClust, iTheta) = mean(spikeCounts{iClust, iTheta}) ;
                sSpikeCount(iClust, iTheta) = std(spikeCounts{iClust, iTheta}) ;
                
                
            end % thetas loop
        end
        %%
        figure2
%        nrows = nClust/4;
        nrows = 3;
        nrows = floor(nrows) + rem(nClust, 4);
        
        for iClust = 1:nClust
            %figure
            %figure(double(clusterID(iClust)+1))
            %subplot(nrows, 4, iClust);
            subplot(nrows, 6, iClust);
            %subplot(nrows, 8, iCh);
            
            hold on
            errorbar(thetaList, mResponseByDir(iClust, :), semResponseByDir(iClust, :), 'o', 'Color', [.7 .7 .7])
            
            %             % TEST TEST just highlighting the directions presented for
            %             % discrimination
            %         errorbar(0, mResponseByDir(iClust, 1), semResponseByDir(iClust, 1), 'o', 'color', [0 0 1])
            %         errorbar(180, mResponseByDir(iClust, 6), semResponseByDir(iClust, 6), 'o', 'color', [1 0 0])
            
            [xvals, respFit] = dsFit(thetaList, mResponseByDir(iClust, :)');
            plot(xvals(1:360),respFit(1:360),'color',[0 0 .6],'linewidth', 2)
            %plot(xvals(1:360),respFit(1:360),'color',[.3 .3 .3],'linewidth', 2)
            
            % save out fits for each cluster, find theta val at max fit.
            allFits(iClust, :)    = respFit;
            [~, maxTheta(iClust)] = max(respFit);
            
            switch dataSource
                case 'kiloSort'
                    %title( num2str(clusterID(iClust)) );
                case 'plx'
                    title(['channel ' num2str(channelNumber)]);
            end
            
            xlabel('Direction')
            ylabel('spks/s')
            xlim([ -10 360])
            %        set(gca, 'XTick', 0:45:359);
            set(gca, 'XTick', 0:90:359);
            %axis square
        end
        % subplot(1,2,2);
        % rho = [mResponseByDir(iClust, :) mResponseByDir(iClust, 1)]; % 'circular' mean response by direction
        % polarplot(thetaListRad, rho, '-o'); hold on
        supertitle(finalDir.name(1:11), 12);
    end %iCh

    %allMaxTheta = [allMaxTheta, maxTheta];
    %clearvars -except baseDir subject sList expts allMaxTheta comp allMaxTheta_end allMaxTheta_start
end % iS


%%
[~, thIx] = sort(maxTheta);

figure
for iClust = 1:nClust
%    subplot(nrows, 4, iClust); hold on
    subplot(nrows, 6, iClust); hold on
    
    errorbar(thetaList, mResponseByDir(thIx(iClust), :), semResponseByDir(thIx(iClust), :), 'o', 'Color', [.7 .7 .7])
    plot(xvals(1:360), allFits(thIx(iClust), 1:360),'color',[0 0 .6],'linewidth', 2)
    
    xlabel('Direction')
    ylabel('spks/s')
    xlim([ -10 360])
    set(gca, 'XTick', 0:90:359);
end
