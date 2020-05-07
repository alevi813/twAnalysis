function [dprime, dirPref] = twagDprime(sumCoh, theseSpikes, centerTime)

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

dprime = (mR - mL)/seDP;
dprime = round(dprime{jj}(ii), 2);
% ^^^ Let's set it up so that:
% positive dprime means right directional preference
% negative dprime means left directional preference
if dprime > 0
    dirPref = 1; %right
else
    dirPref = 0; %left
end

end