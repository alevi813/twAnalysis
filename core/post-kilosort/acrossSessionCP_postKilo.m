fullDir = ('/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/nancy');
sessions = dir(fullDir);
sessions=sessions(~ismember({sessions.name},{'.','..','.DS_Store','allSessions.mat'})); % CLEAN DAT SHIT UP, also keep it from crashing if they already ran it and have a sessions.mat

% sessions = sessions(1:21); % testing

%%
for nSession = 1:length(sessions)
    dateDir = ([fullDir filesep sessions(nSession).name]);
    expts = dir(dateDir);
    expts = expts(~ismember({expts.name},{'.','..','.DS_Store'})); % CLEAN UP
    
    if length(expts) > 1
        expts = expts(1);
    end
    
    %% load  pp
    %  PDS, clock info, & rez
    exDir = ([dateDir filesep expts.name]);
    
    % annoying workaround for loading pp...
    f = dir(exDir);
    f = f(~ismember({f.name},{'.','..','.DS_Store'})); % CLEAN UP
    
    for ii = 1:length(f)
        f(ii).ispp = strcmp(f(ii).name(end-5:end), 'pp.mat');
    end
    ispp = arrayfun(@(x) x.ispp, f, 'uniformoutput', 0);
    ispp = cell2mat(ispp);
    f = f(ispp);
    load([exDir filesep f.name]);
    
    maxScaledKernel = (pp.kernel.b - min(pp.kernel.b))/(max(pp.kernel.b)-min(pp.kernel.b));
    allSession_MaxScaledKernel(:,nSession) = maxScaledKernel;
end

figure(1); subplot(2,1,1)
%imagesc(upsamp(fullKernelBySession,10000)); colorbar; colormap(bone)
imagesc(allSession_MaxScaledKernel); colorbar; colormap(bone)
xlabel('Session')
ylabel('Pulse')
set(gca, 'Ytick', 1:7)
title('Pulse weight')

hold on
x=[10.5, 10.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(x,y,'-r', 'linewidth', 2)

xx=[21.5, 21.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(xx,y,'-g', 'linewidth', 2)

%% CP
path      = '/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData';

% sessionList = [1:21];
sessionList = 1:length(sessions);
subject = 'nancy';


for iSession = 1:length(sessionList)
    
    processedFile = [path filesep 'nancySession' num2str(sessionList(iSession))];
    load([processedFile '.mat'])
    
    allSession_PulseCP(:, iSession) = mPulseCP;
    
    maxScaledPulseCP = (mPulseCP - min(mPulseCP)) / (max(mPulseCP) - min(mPulseCP));
    allSession_maxScaledPulseCP(:, iSession) = maxScaledPulseCP;
end

subplot(2,1,2)
%imagesc(upsamp(fullMaxScaledPulseCP,10000)); colorbar; colormap(bone)
imagesc(allSession_maxScaledPulseCP); colorbar; colormap(bone)
xlabel('Session')
ylabel('Pulse')
set(gca, 'Ytick', 1:7)
title('Pulse CP')

hold on
x=[10.5, 10.5];
y=[0.5 7.5];
plot(x,y,'-r', 'linewidth', 2)

xx=[21.5, 21.5];
%x=[9.5, 9.5];
y=[0.5 7.5];
plot(xx,y,'-g', 'linewidth', 2)