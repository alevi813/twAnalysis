function pp = analyzeBehaviorSession(baseDir, pdsf, trStart, plotOn, ppNum)

%% load PDS, convert from cell, compute PMF & PPK

% baseDir = '/Volumes/LKCLAB/EPHYS/DATA/Projects/nancyTwag2018/PDS';
% pdsf = 'nancy20180914twag.twag1039.PDS';

if nargin < 5
    ppNum = NaN;
    if nargin < 4
        plotOn = false;
        if nargin < 3
            trStart = 1;
        end
    end
end

load([baseDir filesep pdsf], '-mat');

pp = twagConvertFromCell(pdsf, baseDir, trStart);

pp = twagCompute(pp, plotOn);

%save('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/pp/leo-e01-pp.mat', 'pp')
%save('/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/pp/leo-e21-pp.mat', 'pp')

if ~isnan(ppNum)    
    [~, host] = system('hostname');
    
    if strcmp(host(1), 'd') %dhcp-129-116-178-237.cps.utexas.edu
        %desktop
        save(['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/pp/leo-e' num2str(ppNum) '-pp.mat'], 'pp')
    else
        %laptop
        save(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/pp/leo-e' num2str(ppNum) '-pp.mat'], 'pp')
    end    
end

%%

% choFro = nansum(pp.targRchosen(pp.frozenTrials)) / length(pp.targRchosen(pp.frozenTrials));
% disp(choFro);

% 
% 
% for iDist = 1:5
%     
%     dist_sumCoh = pp.sumCoh(pp.stimDistNum==iDist);
%     dist_pulses = pp.cohmat(pp.stimDistNum==iDist, :);
%     
%     figure
%     supertitle(['Dist ' num2str(iDist)], 16);
%     for iPulse = 1:pp.nPulses(1)
%         
%         subplot(pp.nPulses(1), 1, iPulse)
%         scatter(dist_sumCoh, dist_pulses(:,iPulse));
%         lsline;
%         xlabel('full trial net coherence')
%         %xlim([30 70])
%         ylabel(['Pulse ' num2str(iPulse) ' coherence'])
%         %ylim([0 100])
%         
%         r = corrcoef(dist_sumCoh, dist_pulses(:,iPulse));
%         
%         title(['r = ' num2str(r(2))]);
%     end
%     
% end
% 
% 
% %%
% 
% for iDist = 1:5
%     
%     dist_sumCoh = pp.sumCoh(pp.stimDistNum==iDist);
%     dist_pulses = pp.cohmat(pp.stimDistNum==iDist, :);
%     
%     figure
%     for iPulse = 1:pp.nPulses(1)  
%         subplot(pp.nPulses(1), 1, iPulse)
%         histogram(dist_pulses(:,iPulse), 15);
% 
%         xlim([-100 100]);
%     end
% end