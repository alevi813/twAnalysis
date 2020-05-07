function [mSpikesByDist, sSpikesByDist, tSpCnt, dprime, dirPref] = twagPSTH(plotOn, channelNumber, spikeTimes, stimDistNum, centerTime, sumCoh)

% function to calculate and plot PSTH from twag session and calculate d'
% based on firing of R vs L trials (strong, weak, and revco).
%
% INPUTS:
%    - plotOn: 1 if you want to plot, 0 if not
%    - channelNumber: vector of channel numbers from recording that you
%    wish to load.
%    - spikeTimes: nCh x 1 cell with spike times from plx_ts
%    - stimDistNum: from PDS, renamed in twagPlxAnalysis
%    - centerTime: beginning of window for PSTH using trial start and
%    motion onset times, calculated in twagPlxAnalysis
% OUTPUTS:
%    - mSpikesByDist: mean spikes over time during motion by stimulus
%    strength distribution (plotted in PSTH)
%    - sSpikesByDist: std of above
%    - tSpCnt: timecourse of spiking during motion on for individual trials
%    - dprime: d' for all channels
%    - dirPref: direction for which the channel fired more
%       RIGHT = 1, LEFT = 0
% -------------------------------------------------------------------------

% get number of channels and preallocate some vars
nCh           = length(channelNumber);
dprime        = nan(nCh, 1);
dirPref       = nan(nCh, 1);
mSpikesByDist = cell(nCh, 1);
sSpikesByDist = cell(nCh, 1);
tSpCnt        = cell(nCh, 1);

if plotOn == 1
figure; hold on
end

kk = 1;
%for jj = 1:nCh
for jj = channelNumber(kk):channelNumber(nCh)
    % CALCULATE AND PLOT PSTH
    ixix = [1 2 3 4 5];
    clr  = [1 0 0; 1 .6 .6; 0 0 0; .6 .6 1; 0 0 1];
    
    for ii = 1:length(ixix)
        distIx   = stimDistNum(:)==ixix(ii);
        if length(distIx) ~= length(centerTime)     % fix for sessions that run longer after turning of recording
            distIx = distIx(1:length(centerTime));
        end
        [mSpikesByDist{jj}, sSpikesByDist{jj} ,bc, ~, tSpCnt{jj}] = pdsa.eventPsth(spikeTimes{jj, 1}, centerTime(distIx), [-.2 1.5],  .001, ones(100,1)/100);
        
        if plotOn == 1
        %subplot(1, nCh, kk); hold on
        %subplot(6, 4, kk); hold on
        subplot(4, 8, kk); hold on

        plot(bc, mSpikesByDist{jj}, 'Color', clr(ii,:), 'LineWidth', 2)
        xlim([-.2 1.5]);
        %text(-0.4, max(mSpikesByDist{jj})-5, ['ch ' num2str(channelNumber(kk))]);
        end
    end
    
    % CALCULATE d'
    ixLeft  = sumCoh < 0;
    ixRight = sumCoh > 0;
    if length(ixLeft) ~= length(centerTime)     % fix for sessions that run longer after turning of recording
        ixLeft  = ixLeft(1:length(centerTime));
        ixRight = ixRight(1:length(centerTime));
    end
    
    [mLeftTrials,sLeftTrials,~, ~, ~]=pdsa.eventPsth(spikeTimes{jj, 1}, centerTime(ixLeft), [-.2 1.5], .001, ones(100,1)/100);
    [mRightTrials,sRightTrials, ~, ~, ~]=pdsa.eventPsth(spikeTimes{jj, 1}, centerTime(ixRight), [-.2 1.5], .001, ones(100,1)/100);
    
    mL = mean(mLeftTrials); %strong left
    sL = mean(sLeftTrials);
    mR = mean(mRightTrials); %strong right
    sR = mean(sRightTrials);
    
    seDP = sqrt((sL^2 + sR^2)/2);
    
    dprime(jj) = (mR - mL)/seDP;
    dprime(jj) = round(dprime(jj), 2);
    % ^^^ Let's set it up so that:
    % positive dprime means right directional preference
    % negative dprime means left directional preference
    if dprime(jj) > 0
        dirPref(jj) = 1; %right
    else
        dirPref(jj) = 0; %left
    end
    kk = kk+1;
end




