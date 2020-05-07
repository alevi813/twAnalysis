load('/Users/Aaron/twag/dataProcessed/leo/bAll.mat');
load('/Users/Aaron/twag/dataProcessed/leo/summary.mat');
load('/Users/Aaron/twag/dataProcessed/leo/bSession.mat')

exname = arrayfun(@(x) x.exname, bSession, 'UniformOutput', false);

for iCond = 1:3
    figure2
    supertitle(bAll(iCond).temporalWeighting, 13)
    
    [minVal, minIx] = min(summary(iCond).wSlope.dist);
    [maxVal, maxIx] = max(summary(iCond).wSlope.dist);
    
    medVal = median(summary(iCond).wSlope.dist);
    %medIx  = find(summary(iCond).wSlope.dist==medVal);
    [~, medIx] = min(abs(summary(iCond).wSlope.dist-median(summary(iCond).wSlope.dist)));
    
    s_minIx = find(strcmp(exname, bAll(iCond).exname{minIx}));
    s_maxIx = find(strcmp(exname, bAll(iCond).exname{maxIx}));
    s_medIx = find(strcmp(exname, bAll(iCond).exname{medIx}));
    
    % min    
    subplot(1,3,1)
    plot(bSession(s_minIx).ppk)
    title(['min = ' num2str(minVal)])
    axis square

    % med
    subplot(1,3,2)
    plot(bSession(s_medIx).ppk)
    title(['median = ' num2str(medVal)])
    axis square

    % max
    subplot(1,3,3)
    plot(bSession(s_maxIx).ppk)
    title(['max = ' num2str(maxVal)])
    axis square
end

% leoExps = dir('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo/stim');
% leoExps = leoExps(3:end);



