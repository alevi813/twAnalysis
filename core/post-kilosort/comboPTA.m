function comboPTA(subject, condition, spNum)
% get avg PTA from all neurons from a given monkey and stimulus condition 

%%

% determine which sessions to load, given subject and condition
switch subject
    case 'nancy'
        switch condition
            case 'flat'
                sessionList = 1:10; 
            case 'late'
                sessionList = 11:21; 
            case'early'
                sessionList = 22:36; 
        end
    case 'leo'
        switch condition
            case 'flat'
                sessionList = 1:13; 
            case 'late'
                sessionList = 14:24; 
            case'early'
                sessionList = 25:35; 
        end        
end

comp = getComp;

if strcmp(comp, 'laptop')
    loadPath      = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/PTA/';
else
    loadPath      = '/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/PTA/';
end

% initialize a few fields
numSessions = length(sessionList);
combo.PTA  = [];
 
% load session data, concatenate
for ii = 1:numSessions
    processedFile = [loadPath subject 'Session' num2str(sessionList(ii))];
    load([processedFile '.mat'])    
   
    combo.PTA = [combo.PTA, sessionPTA];        
end


% Reshape, take the mean.
combo.PTA  = cat(3, combo.PTA{:});
bigMeanPTA = mean(combo.PTA, 3);

% plot
subplot(3,3,spNum); hold on
clr = gray(8);
%clr = hot(12);
%clr = winter(7);
plot(1:50, bigMeanPTA(:,1), 'Color',  clr(1,:))
plot(51:100, bigMeanPTA(:,2), 'Color',clr(2,:))
plot(101:150, bigMeanPTA(:,3),'Color',clr(3,:))
plot(151:200, bigMeanPTA(:,4),'Color',clr(4,:))
plot(201:250, bigMeanPTA(:,5),'Color',clr(5,:))
plot(251:300, bigMeanPTA(:,6),'Color',clr(6,:))
plot(301:350, bigMeanPTA(:,7),'Color',clr(7,:))


%% plot

% % subplot(2,3,3);
% plot(1:7, maxBigMeanPTA,...
%     'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
% title('Max PTA value by pulse')
% xlabel('Pulse')
% ylabel('Pta max value')
% set(gca, 'Xtick', 1:7)
% xlim([.5 7.5]);

