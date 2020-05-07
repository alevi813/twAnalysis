clear all

subject = 'nancy';

switch subject
    case 'nancy'
        processedDataPath = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData/nancySession';
        
        sessionList = 1:10; tstring = 'baseline'; % baseline all
        %sessionList = 11:21; tstring = 'late'; % late all
        % sessionList = 22:36; tstring = 'early all'; % all early
        
        %sessionList = 1:37; % all
        
    case 'leo'
        processedDataPath = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData/leoSession';
        
        %sessionList = [1:11 13]; tstring = 'baseline';
        sessionList = 14:20; tstring = 'late';
end

conditionCPm = nan(length(sessionList), 2000);
conditionCPs = nan(length(sessionList), 2000);
nUnits = [];
iloop = 1;
for iSession = sessionList(1):sessionList(end)
    load([processedDataPath num2str(iSession) '.mat']);
    
    nUnits = [nUnits goodClustIx];
    
    if ~isempty(cp.allrevco.session.m)
        %     conditionCPm(iSession, :) = mSessionCP_allRevco;
        %     conditionCPs(iSession, :) = sSessionCP_allRevco;
        %     conditionCPm(iloop, :) = cp.allrevco.session.m;
        %     conditionCPs(iloop, :) = cp.allrevco.session.s;
        conditionCPm(iloop, :) = cp.frozenonly.session.m;
        conditionCPs(iloop, :) = cp.frozenonly.session.s;
        
    else
        
    end
    iloop = iloop+1;
end

nUnits = sum(nUnits);

baseline.cpm = nanmean(conditionCPm);
baseline.cps = nanmean(conditionCPs);
baseline.cpsem = baseline.cps/sqrt(length(sessionList));
%baseline.cpsem = baseline.cps/sqrt(nUnits);

figure; hold on
plot(bc, baseline.cpm,'LineWidth', 1.75, 'Color', [0 0 0])
plot(bc, (baseline.cpm + baseline.cpsem), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);
plot(bc, (baseline.cpm - baseline.cpsem), 'LineWidth', 1.75, 'LineStyle', '--', 'Color', [0 0 0]);

title([tstring ' --- nUnits = ' num2str(nUnits)])

h = refline(0, 0.5);
set(h, 'color', 'r')
axis square
ylabel('CP', 'fontsize', 12)
xlabel('time (s)', 'fontsize', 12)
%xlim([-.02 1.2])