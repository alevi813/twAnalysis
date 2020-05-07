function [dprime, dirPref] = twag_dprime(theseSpikes, sumCoh, centerTime, condition)

% CALCULATE d' (let's make this it's own function too!!!)
ixLeft  = sumCoh < 0;
ixRight = sumCoh > 0;
if length(ixLeft) ~= length(centerTime)     % fix for sessions that run longer after turning of recording
    ixLeft  = ixLeft(1:length(centerTime));
    ixRight = ixRight(1:length(centerTime));
end

if sum(theseSpikes) > 0 % whhhhhyyyyyyyyY??
    [mLeftTrials,sLeftTrials,~, ~, ~]=pdsa.eventPsth(theseSpikes, centerTime(ixLeft), [-.5 1.5], .001, ones(100,1)/100);
    [mRightTrials,sRightTrials, ~, ~, ~]=pdsa.eventPsth(theseSpikes, centerTime(ixRight), [-.5 1.5], .001, ones(100,1)/100);
end

if strcmp(condition, 'baseline');
    mL = mean(mLeftTrials(1, 500:1600)); %strong left
    sL = mean(sLeftTrials(1, 500:1600));
    mR = mean(mRightTrials(1, 500:1600)); %strong right
    sR = mean(sRightTrials(1, 500:1600));
elseif strcmp(condition, 'late');
    mL = mean(mLeftTrials(1, 1100:1600)); %strong left
    sL = mean(sLeftTrials(1, 1100:1600));
    mR = mean(mRightTrials(1, 1100:1600)); %strong right
    sR = mean(sRightTrials(1, 1100:1600));
else strcmp(condition, 'early');
    mL = mean(mLeftTrials(1, 500:1000)); %strong left
    sL = mean(sLeftTrials(1, 500:1000));
    mR = mean(mRightTrials(1, 500:1000)); %strong right
    sR = mean(sRightTrials(1, 500:1000));
end
seDP = sqrt((sL^2 + sR^2)/2);

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