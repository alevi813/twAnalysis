%% sliding window, calculating kernel over all trials
% clear all;
% close all;

%load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-allEphys-10162016.mat');

winSize = 500;
nTrials = length(combo.cohmat);
cohmat      = combo.cohmat;
stimDistNum = combo.stimDistNum;
targRchosen = combo.targRchosen;

for ii = 1:nTrials-winSize
    tempPulses      = cohmat(ii:(ii+winSize-1),:);
    tempStimDist    = stimDistNum(ii:(ii+winSize-1));
    tempTargRchosen = targRchosen(ii:(ii+winSize-1));
        
    ppk = ppkTools(tempPulses, tempTargRchosen', 'ridge', true);
    
    fullKernel(:,ii) = ppk.w_norm;
    kernelSlope(ii) = table2array(ppk.w_lm.Coefficients(2,1));

    %some maintenence for finding the transition period switching from bl
    %to late
    if (ii+winSize-1) == 13110
        transitionStartTrial(1) = ii;
    end
    
    if (ii+winSize-1) == 29194
        transitionStartTrial(2) = ii;
    end
end

figure
plot(1:length(kernelSlope), kernelSlope, 'o')

% figure
% subplot(2,2,1)
% imagesc(fullKernel); colorbar; colormap(jet)
% xlabel('Trial start')
% ylabel('Pulse')
% set(gca, 'Ytick', 1:7)
% title('Pulse weight')

% hold on
% x=[12111.5, 12111.5];
% y=[0.5 7.5];
% plot(x,y, 'color', [.5 .5 .5], 'linewidth', 1.5)
% 
% x=[13110.5, 13110.5];
% y=[0.5 7.5];
% plot(x,y,'-k', 'linewidth', 1.5)

%% kernel by session
subject = 'nancy';
%eList   = 57:78;
eList   = [58 60 61 62 64 65 66 67 68 69 70 71 72 74 75 76 77 78 79]; %list of behavior ppfiles that have corresponding ephys
%eList   = [60 61 62 64 65 66 67 68 69 70 71 72 74 75 76 77 78 79];
[ppAll] = ppLoader(subject, eList);

for ii = 1:length(eList)
    fullKernelBySession(:,ii) = ppAll{ii}.pp.kernel.b;
end

subplot(2,2,2)
imagesc(fullKernelBySession); colorbar; colormap(jet)
xlabel('Session')
ylabel('Pulse')
set(gca, 'Ytick', 1:7)
title('Pulse weight')

hold on
x=[10.5, 10.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(x,y,'-k', 'linewidth', 1.5)

%% CP by pulse by session

path      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/';
sessionList = [3 5 6 7 9 10 11 12 13 14 15 16 17 19 20 21 22 23 24];
%sessionList = [5 6 7 9 10 11 12 13 14 15 16 17 19 20 21 22 23 24];
numFlatSessions = length(sessionList);
subject = 'nancy';


for count = 1:numFlatSessions
    processedFile = [path subject 'Session' num2str(sessionList(count))];
    load([processedFile '.mat'])
    
    % raw sesion pulse cp
    allsessionPulseCPm_allRevco(:,count) = sessionPulseCPm_allRevco;
    
    % max scaled pulse cp
    maxScaledPulseCP = (sessionPulseCPm_allRevco - min(sessionPulseCPm_allRevco))/(max(sessionPulseCPm_allRevco)-min(sessionPulseCPm_allRevco));
    fullMaxScaledPulseCP(:,count) = maxScaledPulseCP;
    
    % mean subtracted session pulse cp
    normPulseCP = sessionPulseCPm_allRevco - mean(sessionPulseCPm_allRevco);
    allSessionNormPulseCPm(:,count) = normPulseCP;
    
    % norm normalized session pulse cp
    l2normPulseCP(:,count) = sessionPulseCPm_allRevco/norm(sessionPulseCPm_allRevco);
end

subplot(2,2,4)
imagesc(allSessionNormPulseCPm); colorbar; colormap(jet)
xlabel('Session')
ylabel('Pulse')
set(gca, 'Ytick', 1:7)
title('Pulse CP')

hold on
x=[10.5, 10.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(x,y,'-k', 'linewidth', 1.5)

%%%%%%%%
subplot(2,2,3)
imagesc(fullMaxScaledPulseCP); colorbar; colormap(jet)
xlabel('Session')
ylabel('Pulse')
set(gca, 'Ytick', 1:7)
title('Pulse CP')

hold on
x=[10.5, 10.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(x,y,'-k', 'linewidth', 1.5)


%% CP using sliding window

% channelNumber = 1:24;
% path      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/';
% sessionList = [5 6 7 9 10 11 12 13 14 15 16 17 20 21 22];
% numFlatSessions = length(sessionList);
% subject = 'nancy';
%
%
% for count = 1:numFlatSessions
%     processedFile = [path subject 'Session' num2str(sessionList(count))];
%     load([processedFile '.mat'])
%
%
% end
%


