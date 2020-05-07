%% kernel by session
subject = 'nancy';
%eList   = 57:78;
%eList   = [58 60 61 62 64 65 66 67 68 69 70 71 72 74 75 76 77 78 79 80]; %list of behavior ppfiles that have corresponding ephys
%eList   = [60 61 62 64 65 66 67 68 69 70 71 72 74 75 76 77 78 79 80];
eList   = [58 60 61 62 64 65 66 67 68 69 70 71 72 74 75 76 77 78 79 80 81 82 83 85 86 87 88 89 90 91 92 93 99 100 101]; %list of behavior ppfiles that have corresponding ephys
%eList   = [58 60 61 62 64 65 66 67 68 69 70 71 72  74 75 76 77 78 79 80]; %list of behavior ppfiles that have corresponding ephys

eList   = [104 105 106 107 109 110];

[ppAll] = ppLoader(subject, eList);

for ii = 1:length(eList)
    %     fullKernelBySession(:,ii) = ppAll{ii}.pp.kernel.b;
    %
    %     maxScaledKernel = (ppAll{ii}.pp.kernel.b - min(ppAll{ii}.pp.kernel.b))/(max(ppAll{ii}.pp.kernel.b)-min(ppAll{ii}.pp.kernel.b));
    %     fullMaxScaledKernel(:,ii) = maxScaledKernel;
    fullKernelBySession(:,ii) = ppAll{ii}.pp.ppk.w_norm;
    
    maxScaledKernel = (ppAll{ii}.pp.ppk.w_norm - min(ppAll{ii}.pp.ppk.w_norm))/(max(ppAll{ii}.pp.ppk.w_norm)-min(ppAll{ii}.pp.ppk.w_norm));
    fullMaxScaledKernel(:,ii) = maxScaledKernel;
end

figure(10); subplot(3,1,1)
imagesc(upsamp(fullKernelBySession,10000)); colorbar; colormap(bone)
%imagesc(fullKernelBySession); colorbar; colormap(bone)
xlabel('Session')
ylabel('Pulse')
set(gca, 'Ytick', 1:7)
title('Pulse weight')

hold on
x=[10.5, 10.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(x,y,'-k', 'linewidth', 1.5)

hold on
x=[20.5, 20.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(x,y,'-k', 'linewidth', 1.5)

%%%%%%%%
figure(20); subplot(3,1,1)
imagesc(upsamp(fullMaxScaledKernel,10000)); colorbar; colormap(bone)
%imagesc(fullMaxScaledKernel); colorbar; colormap(bone)
xlabel('Session')
ylabel('Pulse')
set(gca, 'Ytick', 1:7)
title('Pulse weight')

hold on
x=[10.5, 10.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(x,y,'-k', 'linewidth', 1.5)

hold on
x=[20.5, 20.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(x,y,'-k', 'linewidth', 1.5)

%% CP by pulse by session

path      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/';
%sessionList = [3 5 6 7 9 10 11 12 13 14 15 16 17 19 20 21 22 23 24 25];
%sessionList = [5 6 7 9 10 11 12 13 14 15 16 17 19 20 21 22 23 24 25];
sessionList = [3 5 6 7 9 10 11 12 13 14 15 16 17 19 20 21 22 23 24 25 26 27 28 30 31 32 33 34 35 36 37 38 39 40 41];
%sessionList = [3 5 6 7 9 10 11 12 13 14 15 16 17 19 20 21 22 23 24 25];

% path      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/twag/';
% sessionList = [1:21];

path      = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/nForm/';
sessionList = 1:6;

numFlatSessions = length(sessionList);
subject = 'nancy';


for count = 1:numFlatSessions
    processedFile = [path subject 'Session' num2str(sessionList(count))];
    load([processedFile '.mat'])
    
%     % raw sesion pulse cp
%     if count < 3
%         allsessionPulseCPm_allRevco(:,count) = cp.frozenonly.pulse.m(1,:);
%         % max scaled pulse cp
%         maxScaledPulseCP = (cp.frozenonly.pulse.m(1,:) - min(cp.frozenonly.pulse.m(1,:)))/(max(cp.frozenonly.pulse.m(1,:))-min(cp.frozenonly.pulse.m(1,:)));
%     else
%         allsessionPulseCPm_allRevco(:,count) = cp.frozenonly.pulse.m(:,1);
%         % max scaled pulse cp
%         maxScaledPulseCP = (cp.frozenonly.pulse.m(:,1) - min(cp.frozenonly.pulse.m(:,1)))/(max(cp.frozenonly.pulse.m(:,1))-min(cp.frozenonly.pulse.m(:,1)));
%     end
%     
%     fullMaxScaledPulseCP(:,count) = maxScaledPulseCP(:);
    
    mResponsePreferred(count,:) = m(1, :);
    mResponseNull(count,:)      = m(5, :);

    % pre nForm below
    %     % raw sesion pulse cp
    %     allsessionPulseCPm_allRevco(:,count) = sessionPulseCPm_allRevco;
    %
    %     % max scaled pulse cp
    %     maxScaledPulseCP = (sessionPulseCPm_allRevco - min(sessionPulseCPm_allRevco))/(max(sessionPulseCPm_allRevco)-min(sessionPulseCPm_allRevco));
    %     fullMaxScaledPulseCP(:,count) = maxScaledPulseCP;
    %
    %     % mean subtracted session pulse cp
    %     normPulseCP = sessionPulseCPm_allRevco - mean(sessionPulseCPm_allRevco);
    %     allSessionNormPulseCPm(:,count) = normPulseCP;
    %
    %     % norm normalized session pulse cp
    %     l2normPulseCP(:,count) = sessionPulseCPm_allRevco/norm(sessionPulseCPm_allRevco);
    %
    %
    %     %%%% PTA
    %     if count == 1
    %         scaledPTA = (maxPTA - min(maxPTA))/(max(maxPTA) - min(maxPTA));
    %         fullScaledPTA(:,count) = scaledPTA;
    %     else
    %         scaledPTA = (maxPTAm - min(maxPTAm))/(max(maxPTAm) - min(maxPTAm));
    %         fullScaledPTA(:,count) = scaledPTA;
    %     end
end

% figure(10); subplot(3,1,2)
% imagesc(upsamp(allSessionNormPulseCPm,10000)); colorbar; colormap(bone)
% %imagesc(allSessionNormPulseCPm); colorbar; colormap(bone)
% xlabel('Session')
% ylabel('Pulse')
% set(gca, 'Ytick', 1:7)
% title('Pulse CP')
%
% hold on
% x=[10.5, 10.5];
% %x=[9.5, 9.5];
% y=[0.5 7.5];
% plot(x,y,'-k', 'linewidth', 1.5)
%
% hold on
% x=[20.5, 20.5];
% %x=[9.5, 9.5];
% y=[0.5 7.5];
% plot(x,y,'-k', 'linewidth', 1.5)

%%%%%%%%
% figure(20); subplot(3,1,2)
% imagesc(upsamp(fullMaxScaledPulseCP,10000)); colorbar; colormap(bone)
% %imagesc(fullMaxScaledPulseCP); colorbar; colormap(bone)
% xlabel('Session')
% ylabel('Pulse')
% set(gca, 'Ytick', 1:7)
% title('Pulse CP')
% 
% hold on
% x=[10.5, 10.5];
% %x=[9.5, 9.5];
% y=[0.5 7.5];
% plot(x,y,'-k', 'linewidth', 1.5)
% 
% hold on
% x=[20.5, 20.5];
% %x=[9.5, 9.5];
% y=[0.5 7.5];
% plot(x,y,'-k', 'linewidth', 1.5)

%%
figure; subplot(1,2,1)
imagesc(mResponsePreferred)
colorbar
colormap(jet)

subplot(1,2,2)
imagesc(mResponseNull)
colorbar
colormap(jet)
%%
%%%% plot PTA
% subplot(3,1,3)
% imagesc(fullScaledPTA); colorbar; colormap(bone)
% xlabel('Session')
% ylabel('Pulse')
% set(gca, 'Ytick', 1:7)
% title('PTA Peaks')
%
% hold on
% x=[10.5, 10.5];
% %x=[9.5, 9.5];
% y=[0.5 7.5];
% plot(x,y,'-k', 'linewidth', 1.5)
%
% hold on
% x=[20.5, 20.5];
% %x=[9.5, 9.5];
% y=[0.5 7.5];
% plot(x,y,'-k', 'linewidth', 1.5)




%% some scatter plots

% weightIndex = (sum(fullMaxScaledKernel(1:3, :)) - sum(fullMaxScaledKernel(5:7, :))) ./ (sum(fullMaxScaledKernel(1:3, :)) + sum(fullMaxScaledKernel(5:7, :)));
%
% cpIndex = (sum(fullMaxScaledPulseCP(1:3, :)) - sum(fullMaxScaledPulseCP(5:7, :))) ./ (sum(fullMaxScaledPulseCP(1:3, :)) + sum(fullMaxScaledPulseCP(5:7, :)));
%
% figure; subplot(1,2,1); clr = bone(28);
% for ii = 1:length(cpIndex)
%     scatter(cpIndex(ii), weightIndex(ii), 50, clr(ii,:))
%     hold on
% end
% %scatter(cpIndex, weightIndex, 3, clr(ii))
% xlim([-1 1])
% xlabel('CP Index')
% ylim([-1 1])
% ylabel('Weight Index')
%
% subplot(1,2,2)
% ixDiff = weightIndex - cpIndex;
% plot(1:length(ixDiff), ixDiff);
% title('weight index - cp index')

% %%
% % gg for all early
% %load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-Flat-allSessions.mat');
%
% % gg for last six early
% load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-Flat-lastSix.mat');
%
% figure; hold on
% errorbar(1:7, gg.kernel.b, gg.kernel.s.se, ...
%     'o-', 'MarkerFaceColor', [.7 .7 .7], 'linewidth', 2, 'color', [.7 .7 .7]);
%
% % % gg for first four late stim sessions (transition)
% % load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-Late-firstFour.mat');
% %
% % errorbar(1:7, gg.kernel.b, gg.kernel.s.se, ...
% %     'o-', 'MarkerFaceColor', [.3 .3 .3], 'linewidth', 2, 'color', [.3 .3 .3]);
%
% % gg for last six late stim sessions (transition)
% load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/gg/nancy-gg-Late-lastSix.mat');
%
% errorbar(1:7, gg.kernel.b, gg.kernel.s.se, ...
%     'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
%
% %legend('Flat', 'Transition', 'Late');
% legend('Flat', 'Late');
%
% % make it pretty
% xlabel('Pulse', 'FontSize', 16)
% set(gca, 'FontSize', 14)
% ylabel('Weight (normalized)', 'FontSize', 16)
% set(gca, 'Xtick', 1:7)
% set(gca, 'FontSize', 14)
% xlim([.5 7.5]);
% set(gcf, 'color', 'w');
%
% %%
% figure; hold on
% %load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/nancyAllFlat.mat')
% load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/nancyLastSixFlate.mat')
%
% errorbar(1:7, pulseCPm_rc, pulseCPse_rc, ...
%     'o-', 'MarkerFaceColor', [.7 .7 .7], 'linewidth', 2, 'color', [.7 .7 .7]);
%
%
% % load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/nancyFirstFourLate.mat')
% %
% % errorbar(1:7, pulseCPm_rc, pulseCPse_rc, ...
% %     'o-', 'MarkerFaceColor', [.3 .3 .3], 'linewidth', 2, 'color', [.3 .3 .3]);
%
% load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/nancyLastSixLate.mat')
%
% errorbar(1:7, pulseCPm_rc, pulseCPse_rc, ...
%     'o-', 'MarkerFaceColor', [0 0 0], 'linewidth', 2, 'color', [0 0 0]);
%
% title('Choice probability by Pulse')
% h = refline(0,0.5);
% set(h, 'Color', [0 0 0], 'LineStyle', '--');
% xlabel('Pulse')
% set(gca, 'FontSize', 14)
% ylabel('CP')
% set(gca, 'FontSize', 14)
% set(gca, 'Xtick', 1:7)
% xlim([.5 7.5]);
% set(gcf, 'color', 'w');
%
% %legend('Flat', 'Transition', 'Late');
% legend('Flat',  'Late');