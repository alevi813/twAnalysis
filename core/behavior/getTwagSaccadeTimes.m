
% function sacTime = getTwagSaccadeTimes(PDS, screenLatency)

%%%% test sessions
% load('/Users/aaronlevi/Desktop/leoPDS/leo20190806twag.twag1234.PDS', '-mat')
% load('/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160726twag.twag1540.PDS', '-mat')

%%

%%%%%%%% data list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %    dirDataPDS   = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/';
dirDataPDS   = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/';

%%%%%%%%%%% NANCY 
%     pdslist = {'nancy/nancy20160726twag.twag1540.PDS',... 
%         'nancy/nancy20160802twag.twag1356.PDS',...        
%         'nancy/nancy20160803twag.twag1233.PDS',...
%         'nancy/nancy20160804twag.twag1441.PDS',...
%         'nancy/nancy20160812twag.twag1151.PDS',...
%         'nancy/nancy20160823twag.twag1405.PDS',...
%         'nancy/nancy20160825twag.twag1324.PDS',...
%         'nancy/nancy20160831twag.twag1330.PDS',...
%         'nancy/nancy20160906twag.twag1421.PDS',...
%         'nancy/nancy20160908twag.twag1248.PDS',...
%         'nancy/nancy20160909twag.twag1247.PDS',...
%         'nancy/nancy20160913twag.twag1452.PDS',...
%         'nancy/nancy20160915twag.twag1328.PDS',...
%         'nancy/nancy20160920twag.twag1310.PDS',...
%         'nancy/nancy20160921twag.twag1449.PDS',...
%         'nancy/nancy20160923twag.twag1314.PDS',...
%         'nancy/nancy20160928twag.twag1401.PDS',...
%         'nancy/nancy20161005twag.twag1454.PDS',...
%         'nancy/nancy20161013twag.twag1255.PDS',...
%         'nancy/nancy20161014twag.twag1359.PDS',...
%         'nancy/nancy20161018twag.twag1330.PDS',...
%         'nancy/nancy20161019twag.twag1403.PDS',...
%         'nancy/nancy20161020twag.twag1301.PDS',...
%         'nancy/nancy20161021twag.twag1349.PDS',...
%         'nancy/nancy20161027twag.twag1343.PDS',...
%         'nancy/nancy20161028twag.twag1315.PDS',...
%         'nancy/nancy20161101twag.twag1329.PDS',...
%         'nancy/nancy20161102twag.twag1230.PDS',...
%         'nancy/nancy20161103twag.twag1259.PDS',...
%         'nancy/nancy20161109twag.twag1315.PDS',...
%         'nancy/nancy20161110twag.twag1339.PDS',...
%         'nancy/nancy20161111twag.twag1257.PDS',...
%         'nancy/nancy20161115twag.twag1427.PDS',...
%         'nancy/nancy20161130twag.twag1345.PDS',...
%         'nancy/nancy20161201twag.twag1407.PDS',...
%         'nancy/nancy20161202twag.twag1433.PDS'
%     };
% stimFiles = dir('/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/nancy/stim');
% stimFiles = stimFiles(4:end);

% %%%%%%%%%%%%%% LEO 
pdslist = {'leo/PDS/leo20190806twag.twag1234.PDS',...
    'leo/PDS/leo20190807twag.twag1310.PDS',...
    'leo/PDS/leo20190814twag.twag1206.PDS',...
    'leo/PDS/leo20190815twag.twag1314.PDS',...
    'leo/PDS/leo20190816twag.twag1228.PDS',...
    'leo/PDS/leo20190827twag.twag1159.PDS',...
    'leo/PDS/leo20190829twag.twag1200.PDS',...
    'leo/PDS/leo20190830twag.twag1219.PDS',...
    'leo/PDS/leo20190904twag.twag1147.PDS',...
    'leo/PDS/leo20190905twag.twag1253.PDS',...
    'leo/PDS/leo20190916twag.twag1218.PDS',...
    'leo/PDS/leo20190917twag.twag1223.PDS',...
    'leo/PDS/leo20190919twag.twag1209.PDS',...
    'leo/PDS/leo20191008twag.twag1148.PDS',...
    'leo/PDS/leo20191009twag.twag1203.PDS',...
    'leo/PDS/leo20191010twag.twag1149.PDS',...
    'leo/PDS/leo20191016twag.twag1235.PDS',...
    'leo/PDS/leo20191021twag.twag1310.PDS',...
    'leo/PDS/leo20191023twag.twag1230.PDS',...
    'leo/PDS/leo20191025twag.twag1323.PDS',...
    'leo/PDS/leo20191120twag.twag1225.PDS',...
    'leo/PDS/leo20191122twag.twag1321.PDS',...
    'leo/PDS/leo20191126twag.twag1305.PDS',...
    'leo/PDS/leo20191129twag.twag1348.PDS',...
    'leo/PDS/leo20191205twag.twag1338.PDS',...
    'leo/PDS/leo20191206twag.twag1310.PDS',...
    'leo/PDS/leo20191207twag.twag1237.PDS',...
    'leo/PDS/leo20191209twag.twag1233.PDS',...
    'leo/PDS/leo20191210twag.twag1400.PDS',...
    'leo/PDS/leo20191211twag.twag1320.PDS',...
    'leo/PDS/leo20191212twag.twag1243.PDS',...
    'leo/PDS/leo20191213twag.twag1258.PDS',...
    'leo/PDS/leo20191216twag.twag1309.PDS',...
    'leo/PDS/leo20191217twag.twag1432.PDS',... %%% *TWO* pds files from this day
    'leo/PDS/leo20191218twag.twag1255.PDS';
   
   };
stimFiles = dir('/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo/stim');
stimFiles = stimFiles(4:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nFiles = length(pdslist);
assert (length(stimFiles)==nFiles, 'n stim files does not match n pds files')


%% stimulus timings
figure

%screenLatency = 0.0524634; % thats 52.4624ms
screenLatency = 0;


for iFile = 1:nFiles
%for iFile = 5
    
    load([dirDataPDS pdslist{iFile}], '-mat');
    
    sacTime    = nan(length(PDS.data), 1);
    sacIx    = nan(length(PDS.data), 1);
    fpoffTime_eyelink = nan(length(PDS.data), 1);
    froIx = cell2mat(cellfun(@(x) x.stimulus.frozenTrials, PDS.data, 'UniformOutput', false));
    
    sRate    = PDS.initialParametersMerged.eyelink.srate;
    targHoldDur = sRate*(PDS.initialParametersMerged.stimulus.durationtarghold-.1);

    
    for ii = 1:length(PDS.data)
        if ~isfield(PDS.data{ii}, 'correct')
            PDS.data{ii}.correct = nan;
        end
    end
    correct = cell2mat(cellfun(@(x) x.correct, PDS.data, 'UniformOutput', false));

    for iTrial = 1:length(PDS.data)
        
        if PDS.data{iTrial}.pldaps.goodtrial ~= 1 %|| iTrial == 1121
            % do nothing
            fpoffTime(iTrial)         = NaN;
            sacTime(iTrial)           = NaN;
            sacIx(iTrial)             = NaN;
            fpoffTime_eyelink(iTrial) = NaN;
            
            x_tracePreSac(iTrial, :)  = nan(1,sRate+1);
            y_tracePreSac(iTrial, :)  = nan(1,sRate+1);
        else
            % get fpOffTime from PDS file
            fpoffTime(iTrial)    = PDS.data{iTrial}.stimulus.statesStartTime(6)+ screenLatency;
            trialEndTime(iTrial) = PDS.data{iTrial}.stimulus.statesStartTime(8)+ screenLatency;
            
            %% eyelink data & timings
            % convert eyelink samples to time
            eyetime = (1:PDS.data{iTrial}.eyelink.sampleNum) / sRate;
            
            %%---- eyelink traces
            eyex = PDS.data{iTrial}.eyelink.samples(4,:);
            eyey = PDS.data{iTrial}.eyelink.samples(5,:);
            
            %%---- round to 3pts and make sure divisible by 2... workaround
            fpRound     = round(fpoffTime(iTrial), 3, 'decimals');
            trEndRound     = round(trialEndTime(iTrial), 3, 'decimals');
            
            if mod(fpRound, 1/sRate) > 0
                fpRound = fpRound - .001;
            end
            if mod(trEndRound, 1/sRate) > 0
                trEndRound = trEndRound - .001;
            end
            
            % get the eyelink sample nearest to go signal by checking the
            % rounded fpOff time against the eyelink time
            try
                fpoffEyelinkSample = find(eyetime == fpRound);
                trEndEyelinkSample = find(eyetime == trEndRound);
            catch ME
            end
            
            if isempty(fpoffEyelinkSample)
                [~, fpoffEyelinkSample] = min(abs(eyetime-fpRound));
            end
            assert(~isempty(fpoffEyelinkSample), 'No fpoffsample')
            
            if isempty(trEndEyelinkSample)
                [~, trEndEyelinkSample] = min(abs(eyetime-trEndRound));
            end
            assert(~isempty(trEndEyelinkSample), 'No trial end sample')
            
            % trial time from go signal
            %eyetimeFromGo  = eyetime(fpoffEyelinkSample:end);
            eyetimeFromGo  = eyetime(fpoffEyelinkSample:trEndEyelinkSample);
            fpoffTime_eyelink(iTrial) = eyetime(fpoffEyelinkSample);
            
            % eye trace from go signal, diffed
%             fromGoEyex = eyex(fpoffEyelinkSample:end);            
%             fromGoEyey = eyey(fpoffEyelinkSample:end);
            fromGoEyex = eyex(fpoffEyelinkSample:trEndEyelinkSample);
            fromGoEyey = eyey(fpoffEyelinkSample:trEndEyelinkSample);
            
            %fromGoEyex = diff(fromGoEyex);
            %fromGoEyey = diff(fromGoEyey);
            
            %         trace = [fromGoEyex; fromGoEyey];
            %         result = saccadeDetector(fromGoTime, trace);
            
            %% ---- find saccade
            %[sacTime(iTrial), sacIx(iTrial)] = findSaccade(eyetimeFromGo, [fromGoEyex; fromGoEyey]);
            sacIx(iTrial) = fpWindowExit([fromGoEyex; fromGoEyey], PDS.initialParametersMerged.stimulus.fpRect, sRate);
            if sacIx(iTrial)==Inf
%                 sacTime(iTrial) = NaN;                 
%                 x_tracePreSac(iTrial, :) = nan(size(x_tracePreSac(1,:)));
%                 y_tracePreSac(iTrial, :) = nan(size(x_tracePreSac(1,:)));

                % if you couldn't find a sacIx/Time, just save it as the go signal time
                sacTime(iTrial) = eyetimeFromGo(1); % if you couldn't find a sacIx/Time, just save it as the go signal time
            else
                sacTime(iTrial) = eyetimeFromGo(sacIx(iTrial));
            end    
            
            x_tracePreSac(iTrial, :) = eyex(find(eyetime==sacTime(iTrial))-sRate : find(eyetime==sacTime(iTrial)) );
            y_tracePreSac(iTrial, :) = eyey(find(eyetime==sacTime(iTrial))-sRate : find(eyetime==sacTime(iTrial)) );
            %end
            
            
        end % goodtrial check
    end % trial loop
    
    relativeSacTime = sacTime - fpoffTime_eyelink;
    
%     figure(1)
%     subplot(5, ceil(nFiles/5), iFile)
%     histogram(relativeSacTime)
%     %histogram(relativeSacTime(~froIx)); hold on
%     %histogram(relativeSacTime(froIx))
% %     histogram(relativeSacTime(correct==1)); hold on
% %     histogram(relativeSacTime(correct==0))
% 
%     figure(2) 
%     subplot(5, ceil(nFiles/5), iFile)
%     plot(nanmean(x_tracePreSac)); hold on; yl = ylim;
%     plot([sRate sRate], [yl(1) yl(2)], 'k--');
% 
%     figure(3) 
%     subplot(5, ceil(nFiles/5), iFile)
%     plot(nanmean(y_tracePreSac)); hold on; yl = ylim;
%     plot([sRate sRate], [yl(1) yl(2)], 'k--');
%     
%     figure(4) 
%     subplot(5, ceil(nFiles/5), iFile)
%     plot(nanstd(x_tracePreSac)); hold on; yl = ylim;
%     plot([sRate sRate], [yl(1) yl(2)], 'k--');
%     
%     figure(5) 
%     subplot(5, ceil(nFiles/5), iFile)
%     plot(nanstd(y_tracePreSac)); hold on; yl = ylim;
%     plot([sRate sRate], [yl(1) yl(2)], 'k--');
% 
%     figure(6)
%     plot([PDS.initialParametersMerged.stimulus.fpRect(1) PDS.initialParametersMerged.stimulus.fpRect(3)], [PDS.initialParametersMerged.stimulus.fpRect(4) PDS.initialParametersMerged.stimulus.fpRect(4)], 'k'); hold on
%     plot([PDS.initialParametersMerged.stimulus.fpRect(1) PDS.initialParametersMerged.stimulus.fpRect(3)], [PDS.initialParametersMerged.stimulus.fpRect(2) PDS.initialParametersMerged.stimulus.fpRect(2)], 'k')
%     plot([PDS.initialParametersMerged.stimulus.fpRect(1) PDS.initialParametersMerged.stimulus.fpRect(1)], [PDS.initialParametersMerged.stimulus.fpRect(2) PDS.initialParametersMerged.stimulus.fpRect(4)], 'k')
%     plot([PDS.initialParametersMerged.stimulus.fpRect(3) PDS.initialParametersMerged.stimulus.fpRect(3)], [PDS.initialParametersMerged.stimulus.fpRect(2) PDS.initialParametersMerged.stimulus.fpRect(4)], 'k')
%     plot(x_tracePreSac', y_tracePreSac')
%     xlim([0 PDS.initialParametersMerged.display.pWidth])
%     ylim([0 PDS.initialParametersMerged.display.pHeight])

    
    figure
    plot([PDS.initialParametersMerged.stimulus.fpRect(1) PDS.initialParametersMerged.stimulus.fpRect(3)], [PDS.initialParametersMerged.stimulus.fpRect(4) PDS.initialParametersMerged.stimulus.fpRect(4)], 'k'); hold on
    plot([PDS.initialParametersMerged.stimulus.fpRect(1) PDS.initialParametersMerged.stimulus.fpRect(3)], [PDS.initialParametersMerged.stimulus.fpRect(2) PDS.initialParametersMerged.stimulus.fpRect(2)], 'k')
    plot([PDS.initialParametersMerged.stimulus.fpRect(1) PDS.initialParametersMerged.stimulus.fpRect(1)], [PDS.initialParametersMerged.stimulus.fpRect(2) PDS.initialParametersMerged.stimulus.fpRect(4)], 'k')
    plot([PDS.initialParametersMerged.stimulus.fpRect(3) PDS.initialParametersMerged.stimulus.fpRect(3)], [PDS.initialParametersMerged.stimulus.fpRect(2) PDS.initialParametersMerged.stimulus.fpRect(4)], 'k')
    plot(x_tracePreSac(froIx,:)', y_tracePreSac(froIx)')
    xlim([0 PDS.initialParametersMerged.display.pWidth])
    ylim([0 PDS.initialParametersMerged.display.pHeight])
    
    % append saccade times to stim file and save
    stim         = load([stimFiles(iFile).folder filesep stimFiles(iFile).name]);
    if isfield(stim, 'xTrials')
        stim.sacTime = relativeSacTime(1:stim.xTrials); 
    else
        stim.sacTime = relativeSacTime;
    end
    save([stimFiles(iFile).folder filesep stimFiles(iFile).name], '-struct', 'stim', '-v7.3')  
    
    %sessionSaccades_leo_ajl{iFile} = relativeSacTime; 
%     saccadeTracesX{iFile} = x_tracePreSac;
%     saccadeTracesY{iFile} = y_tracePreSac;
    
    clearvars -except pdslist dirDataPDS nFiles screenLatency stimFiles saccadeTracesX saccadeTracesY %x_tracePreSac y_tracePreSac %sessionSaccades sessionSaccades_ajl sessionSaccades_leo sessionSaccades_leo_ajl
end %iFile