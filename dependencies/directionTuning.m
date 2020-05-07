%% plot and calculate mtrfmap data
% subject and experiment session #

clear all

subject      = 'nancy';   % monkey name
eList        = 5;        % session number. See ephysDataFactory for list.
%eList        = [2 3 4 5 6 7 8 10 12 13 14 16 17 18 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 36 37];
saveSession  = false;
forceProcess = true;

for kk = 1:length(eList)
    clearvars -except subject eList saveSession forceProcess kk
    
    %% check if file already exists for this data set
    savePath      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/mtrfmap/';
    processedFile = [savePath subject 'Session' num2str(eList(kk))];
    
    if exist([processedFile '.mat'], 'file') && forceProcess == false;
        
        disp('This file already exists. Loading...')
        load([processedFile '.mat'])
        disp('Done')
    else
        %%
        [pdsfile, plxfile, PDS, flagPlexonTrialStartTime] = mtrfmapDataFactory(subject, eList(kk));
        
        channelNumber = [1:24];
        nCh = length(channelNumber);
        
        spikeTimes  = cell(nCh, 1);
        goodChannelIx = nan(nCh,1);
        for jj = 1:nCh
            [n, spikeTimes{jj, 1}] = plx_ts(plxfile, channelNumber(jj), 1);% FILENAME, CHANNEL, UNIT
            
            if spikeTimes{jj} > 0
                goodChannelIx(jj) = 1;
            else
                goodChannelIx(jj) = 0;
            end
        end
         goodChannel = find(goodChannelIx);
        
        % PDS clean up... convert necessary vars from cell
        % goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
        % goodtrial = cell2mat(goodtrial);
        % goodtrial(isnan(goodtrial)) = 0;
        
        stimOn = cellfun(@(x) x.stimulus.timeStimOn, PDS.data, 'UniformOutput', false);
        stimOn = cell2mat(stimOn)';
        centerTime = flagPlexonTrialStartTime + stimOn;
        
        nThetas     = 8;
        thetaLength = 3000/nThetas;
        
        thetaChangeTimes = nan(length(PDS.data), nThetas-1);
        for jj = 1:length(PDS.data)
            for ii = 1:nThetas-1
                thetaChangeTimes(jj,ii) = (centerTime(jj)+(thetaLength*ii));
            end
        end
        
        allMotionTimes = [centerTime, thetaChangeTimes];
        
        validIx = ~isnan(allMotionTimes);
        validIx = validIx(:,1);
        
        for ii = 1:length(PDS.data)
            if validIx(ii) == 1
                allThetas(ii,:) = PDS.data{ii}.stimulus.motion1.thetas;
            end
        end
        
        
        
        
        
        %% plot
         figure; 
        thetas    = 0:45:360;
        thetasRad = deg2rad(thetas);
        window    = [0 .375]; %in seconds 0 to .375 is default, is length of each theta presentation
        %window    = [.05 .375];
        
       for jj = 1:length(goodChannel)
%        for jj = 1:24
            for ii = 1:nThetas
                dirIx = allThetas==thetas(ii);
                dirIx = dirIx(:);
                
                [mSpikesByDir{jj}, sSpikesByDir{jj} ,bc, ~, ~] = pdsa.eventPsth(spikeTimes{goodChannel(jj), 1}, allMotionTimes(dirIx), window,  .001, ones(100,1)/100);
%                [mSpikesByDir{jj}, sSpikesByDir{jj} ,bc, ~, ~] = pdsa.eventPsth(spikeTimes{jj}, allMotionTimes(dirIx), window,  .001, ones(100,1)/100);

                %figure; plotRaster(spikeTimes{goodChannel(jj), 1}, allMotionTimes(dirIx), window)
                
                %     subplot(1,2,1); hold on
                %     plot(bc, mSpikesByDir,'LineWidth', 2);
                
                mResponseByDir{jj}(ii) = mean(mSpikesByDir{jj});
            end
            
            rho{jj} = [mResponseByDir{jj} mResponseByDir{jj}(1)]; % 'circular' mean response by direction
                                                                  % i.e. just tacked on the first value on the
                                                                  % end for symmetrical 0-360  
           %subplot(1,length(goodChannel),jj);
           %figure
           %polarplot(thetasRad, rho{jj}, 'r'); hold on
           
           
           % take vector average and plot it in blue on the polar plot
           [x,y] = pol2cart(thetasRad, rho{jj});
           vecAv = [mean(x) mean(y)];
           [cellVectorPref{jj}(1), cellVectorPref{jj}(2)] = cart2pol(vecAv(1), vecAv(2));
           %polarplot(cellVectorPref{jj}(1), cellVectorPref{jj}(2), 'b*')
           
           %subplot(1,length(goodChannel),jj);
           figure
           plot(thetas, rho{jj}, '-o'); hold on
           [xp, fit, s] = fitCos(thetas, rho{jj}); %fit a cosine
           plot(xp,fit(s,xp))
           
           title(['ch ' num2str(goodChannel(jj))]);
        end
        %% save it
        if saveSession == true
            disp('Saving mapping data...');
            save([savePath subject 'Session' num2str(eList(kk))]);
            disp('Session saved');
        end
    end
end