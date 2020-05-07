processedDataDir = '/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData';
sessionList         = dir(processedDataDir);
sessionList         = sessionList(5:end); % hardcoded for now.

listOrder = [1 12 23 32 33 34 35 36 37 2 3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20 21 22 24 25 26 27 28 29 30 31];

sessionList = sessionList(listOrder);

kernelSlope = nan(length(sessionList), 1);
cpSlope     = nan(length(sessionList), 1);

for iSession = 1:length(sessionList)-1 %temporarily leaving out the final session until I combine the two files with plxutil
    %for iSession = 16
    
    session     = load([processedDataDir filesep sessionList(iSession).name]);
    kernelSlope(iSession) = table2array(session.pp.ppk.w_lm.Coefficients(2,1));
    
    if isfield(session.cp.frozenonly.session, 'pulse')
        if ~isnan(session.cp.frozenonly.session.pulse.m)
            cplm    = fitlm(1:7, session.cp.frozenonly.session.pulse.m);
            cpSlope(iSession) = table2array(cplm.Coefficients(2,1));
        else
            cpSlope(iSession) = NaN;
        end
    end
    
    if iSession < 22
        if iSession < 11
            clr = [19/255 40/255 230/255];
        else
            clr = [255/255 179/255 0/255];
        end
    else
        clr = [214/255 4/255 0/255];
        
    end
    
    % plot iteratively
%     figure(1); hold on
%     plot(kernelSlope(iSession), cpSlope(iSession), 'o', 'Color', clr)
    
    figure(2);
    subplot(2,1,1); hold on
    plot(iSession, kernelSlope(iSession), 'o', 'Color', clr)
    
    subplot(2,1,2); hold on
    plot(iSession, cpSlope(iSession), 'o', 'Color', clr)
end


figure(1)
xlabel('kernel slope', 'FontSize', 12);
ylabel('CPt slope', 'FontSize', 12);
xlim([ -.15 .15])
ylim([ -.035 .035])
plot([ -.15 .15], [0 0], 'k--')
plot([0 0], [ -.05 .05], 'k--')

figure(2)
subplot(2,1,1);
ylim([-.145 .145])
ylabel('kernel slope', 'FontSize', 12);
xlim([0 38])
plot([0 38], [0 0], 'k--')

subplot(2,1,2);
ylabel('CPt slope', 'FontSize', 12);
xlabel('Session #', 'FontSize', 12);
xlim([0 38])
ylim([-.035 .035])
plot([0 38], [0 0], 'k--')

%%
figure(4)
subplot(2,1,1)

% plot(1:10,  kernelSlope(1:10), 'Color', [19/255 40/255 230/255]); hold on
% plot(10:21, kernelSlope(10:21), 'Color', [255/255 179/255 0/255])
% plot(21:37, kernelSlope(21:end), 'Color', [214/255 4/255 0/255])


plot(1:10,   kernelSlope(1:10),'-o', 'Color', [19/255 40/255 230/255]); hold on
plot(10:21, kernelSlope(10:21),'-o', 'Color', [255/255 179/255 0/255])
plot(21:37,  kernelSlope(21:end),'-o', 'Color', [214/255 4/255 0/255])

% hold on
% plot(1:10,  cpSlope(1:10), 'Color', [19/255 40/255 230/255], 'LineStyle', '--')
% plot(10:21, cpSlope(10:21), 'Color', [255/255 179/255 0/255], 'LineStyle', '--')
% plot(21:37, cpSlope(21:end), 'Color', [214/255 4/255 0/255], 'LineStyle', '--')

hold on
plot(1:10,   cpSlope(1:10), '--o', 'Color', [19/255 40/255 230/255])
plot(10:21,  cpSlope(10:21), '--o', 'Color', [255/255 179/255 0/255])
plot(21:37,  cpSlope(21:end), '--o', 'Color', [214/255 4/255 0/255])


ylim([-.145 .145])
ylabel('Slope of linear fit', 'FontSize', 12);
xlim([0 38])
%plot([0 38], [0 0], 'k--')





