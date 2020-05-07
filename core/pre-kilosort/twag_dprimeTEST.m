function [dprime, dirPref] = twag_dprimeTEST(theseSpikes, sumCoh, centerTime, condition)

% CALCULATE d' (let's make this it's own function too!!!)
ixLeft  = sumCoh < 0;
ixRight = sumCoh > 0;
if length(ixLeft) ~= length(centerTime)     % fix for sessions that run longer after turning of recording
    ixLeft  = ixLeft(1:length(centerTime));
    ixRight = ixRight(1:length(centerTime));
end

if sum(theseSpikes) > 0 % whhhhhyyyyyyyyY??
    [mLeftTrials,sLeftTrials,~, ~, spkL]=pdsa.eventPsth(theseSpikes, centerTime(ixLeft), [-.5 1.5], .001, ones(100,1)/100);
    [mRightTrials,sRightTrials, ~, ~, spkR]=pdsa.eventPsth(theseSpikes, centerTime(ixRight), [-.5 1.5], .001, ones(100,1)/100);
end

if strcmp(condition, 'baseline')
    mL = nanmean(sum(spkL,2));
    sL = nanvar(sum(spkL,2));
    mR = nanmean(sum(spkR,2));
    sR = nanvar(sum(spkR,2));
elseif strcmp(condition, 'late')
    mL = nanmean(sum(spkL(:, 1100:1600),2));
    sL = nanvar(sum(spkL(:, 1100:1600),2));
    mR = nanmean(sum(spkR(:, 1100:1600),2));
    sR = nanvar(sum(spkR(:, 1100:1600),2));
else strcmp(condition, 'early')
    mL = nanmean(sum(spkL(:, 500:1000),2));
    sL = nanvar(sum(spkL(:, 500:1000),2));
    mR = nanmean(sum(spkR(:, 500:1000),2));
    sR = nanvar(sum(spkR(:, 500:1000),2));
end
seDP = sqrt((sL + sR)/2);

dprime = (mR - mL)/seDP;
dprime = round(dprime, 2);
% ^^^ Let's set it up so that:
% positive dprime means right directional preference
% negative dprime means left directional preference
if dprime > 0
    dirPref = 1; %right
else
    dirPref = 0; %left
end