% pds to stim
% Originally made for nancy flat-stimulus files for Jake & Klaus
% to match Jake's data structure. Updated to use neuroGLM w/ ajl MT data
% under multiple temporal stimulus conditions.
%
% ajl 2018


    %%%%%%%% data list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dirDataPDS   = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/';
    dirDataPDS   = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/';

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% pds 2 stim

for iSession = 26     
    
    load([dirDataPDS pdslist{iSession}], '-mat');
    
    condition = PDS.initialParametersMerged.gabs.condition;
    if strcmp(condition, 'baseline')
        condition = 'flat';
    end
    
    %luminanceMinMuMax = PDS.initialParametersMerged.display.gamma.fullPathLaod.dv.disp.gamma.minMuMax; % wow.
    luminanceMinMuMax = [0.5695 26.5208 52.8500]';  %%% hardcoding rig 1 params % make sure to change for rig 2 %%%
    
    goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
    goodtrial = cell2mat(goodtrial);
    goodtrial(isnan(goodtrial)) = 0;
    goodtrial = logical(goodtrial');
    
    for iTrial = 1:length(PDS.data)
        basephase(iTrial,:,:) = PDS.data{iTrial}.theseGabs.basephase;
        maskphase(iTrial,:,:) = PDS.data{iTrial}.theseGabs.maskphase;
        
        %targchosen
        if ~isfield(PDS.data{iTrial}, 'targRchosen')
            PDS.data{iTrial}.targRchosen = nan;
        end
        
        %correct
        if ~isfield(PDS.data{iTrial}, 'correct')
            PDS.data{iTrial}.correct = nan;
        end
        
        %targcorrect
        if ~isfield(PDS.data{iTrial}.theseGabs, 'targCorrect')
            PDS.data{iTrial}.theseGabs.targCorrect = nan;
        end
        
        % unfortunately, our cohidx is binary (0 not pulsing, 1 pulsing) and we
        % would like it to be signed for if pulse was left (-) or right (+).
        % fix that here
        tmpcohix = double(PDS.data{iTrial}.theseGabs.cohidx);
        
        for iPulse = 1:PDS.data{iTrial}.theseGabs.nPulses
            % if that pulse was to the left (negative), then change cohidx
            if PDS.data{iTrial}.theseGabs.gabCoh(iPulse) < 0
                tmpcohix(iPulse, tmpcohix(iPulse, :) == 1)  =  -1;
            end % is gabCoh negative
            
        end % loop through 7 pulses
        
        pulses(iTrial,:, :) = tmpcohix;
        
    end % loop through trials
    
    nTrials = iTrial;
    nPulses = PDS.data{1}.theseGabs.nPulses;
    nGabors = size(PDS.data{1}.Gpars, 2);
    nFrames = 63;
    
    correct = cellfun(@(x) x.correct, PDS.data, 'UniformOutput', false);
    correct = cell2mat(correct');
    
    targcorrect = cellfun(@(x) x.theseGabs.targCorrect, PDS.data, 'UniformOutput', false);
    targcorrect = cell2mat(targcorrect);
    % right = 2
    targcorrect(targcorrect==1) = 2;
    % left = 1
    targcorrect(targcorrect==0) = 1;
    targcorrect = targcorrect';
    
    
    targchosen  = cellfun(@(x) x.targRchosen, PDS.data, 'UniformOutput', false);
    targchosen  = cell2mat(targchosen);
    % right = 2
    targchosen(targchosen==1) = 2;
    % left = 1
    targchosen(targchosen==0) = 1;
    targchosen = targchosen';
    
    % in the old 'dirprob' for pulse generation, 0 = revco
    % in the stimDistNum version, 3 = revco. Here we conform to dirprob
    dirprob = cellfun(@(x) x.theseGabs.stimDistNum, PDS.data, 'UniformOutput', false);
    dirprob = cell2mat(dirprob');
    dirprob(dirprob==3) = 0;
    dirprob(dirprob==6) = 0;
    
    %    exname = [PDS.initialParametersMerged.session.file(1) PDS.initialParametersMerged.session.file(6:13)];
    exname = [PDS.initialParametersMerged.session.file(1) PDS.initialParametersMerged.session.file(4:11)];
    
    frate = PDS.initialParametersMerged.display.frate;
    
    gaborXY = PDS.data{1}.stimulus.pos';
    
    ppd = PDS.initialParametersMerged.display.ppd;
    
    sf = PDS.initialParametersMerged.stimulus.sf;
    tf = PDS.initialParametersMerged.stimulus.tf;
    
    theta = PDS.initialParametersMerged.stimulus.theta;
    
    targ1XY = PDS.initialParametersMerged.stimulus.targL ./ ppd;
    targ2XY = PDS.initialParametersMerged.stimulus.targR ./ ppd;
    
    targ1XY = repmat(targ1XY, nTrials, 1);
    targ2XY = repmat(targ2XY, nTrials, 1);
    
    trialnumber = cell2mat(cellfun(@(x) x.trialnumber, PDS.data, 'UniformOutput', false));
    
    eyepos = cellfun(@(x) x.eyelink.samples, PDS.data, 'UniformOutput', false);
    
    % trialId
    %%%%%%%%%%%%%%%%%%%%%%
    for iTrial = 1:nTrials
        for jTrial = 1:nTrials
            trialMatrix(iTrial, jTrial) = double(isequal(PDS.data{iTrial}.theseGabs.basephase, PDS.data{jTrial}.theseGabs.basephase));
        end
        
        trialMatrix(iTrial, :) = trialMatrix(iTrial, :) .* iTrial;
        trIx = find(trialMatrix(iTrial, :) > 0);
        
        trialId(trIx) = trialMatrix(iTrial, trIx);
    end
    
    % change the trial id to be the first instance of the frozen trial, instead
    % of the last. This is necessary because of the clunky loop.
    trialMode = mode(trialId);
    trialId(trialId == trialMode) = find(trialId == trialMode, 1);
    
    % frozenTrialIds
    %%%%%%%%%%%%%%%%%%%%%%
    doPlot = false;
    
    if doPlot
        figure;
        subplot(2,1,1); plot(trialId); title(PDS.initialParametersMerged.session.file)
        subplot(2,1,2); histogram(trialId, 66666);
    end
    
    % all these sessions should only have one frozen trial -- so taking the
    % mode should be fine. check the plots above to make sure you don't have
    % more than one frozenId
    frozenTrialIds = mode(trialId);
    
    
    % trialCnt
    %%%%%%%%%%%%%%%%%%%%%%
    %    for iTrial = 1:length(unique(trialId))
    for iTrial = 1:length(trialId)
        trialCnt(iTrial) = sum(trialId == iTrial);
    end
    
    % timing struct
    %%%%%%%%%%%%%%%%%%%%%%
    %load(['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/processedEphys/twag/clock/nancySession' num2str(iSession)])
    %screenLatency = 0.0524634; % thats 52.4624ms
    
%    load(['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo/clock/leoSessionClock' num2str(iSession) '.mat'])
    load(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo/clock/leoSessionClock' num2str(iSession) '.mat'])
    screenLatency = 0;
    
    if length(syncClock.plexonTrialStartTime) ~= nTrials
        xTrials = length(syncClock.plexonTrialStartTime);
    else
        xTrials = nTrials;
    end
    
    %    sacTime = getTwagSaccadeTimes(PDS, screenLatency);
    
    for iTrial = 1:xTrials
        timing(iTrial).fpon        = PDS.data{iTrial}.stimulus.statesStartTime(1)+ screenLatency;
        timing(iTrial).fpentered   = PDS.data{iTrial}.stimulus.statesStartTime(2)+ screenLatency;
        timing(iTrial).targson     = PDS.data{iTrial}.stimulus.statesStartTime(3)+ screenLatency;
        timing(iTrial).motionon    = PDS.data{iTrial}.stimulus.statesStartTime(4)+ screenLatency;
        timing(iTrial).motionoff   = PDS.data{iTrial}.stimulus.statesStartTime(5)+ screenLatency;
        timing(iTrial).fpoff       = PDS.data{iTrial}.stimulus.statesStartTime(6)+ screenLatency;
        %timing(iTrial).saccade     = sacTime(iTrial);
        timing(iTrial).reward      = PDS.data{iTrial}.stimulus.statesStartTime(8)+ screenLatency;
        timing(iTrial).targsoff    = PDS.data{iTrial}.stimulus.statesStartTime(8)+ screenLatency;
        timing(iTrial).breakfix    = PDS.data{iTrial}.stimulus.statesStartTime(7)+ screenLatency;
        
        timing(iTrial).duration = PDS.data{iTrial}.trialend;
        
        plength = (PDS.data{iTrial}.stimulus.statesStartTime(5) - PDS.data{iTrial}.stimulus.statesStartTime(4)) / 7 ;
        timing(iTrial).pulses(1) = timing(iTrial).motionon;
        for ipulse = 1:6
            timing(iTrial).pulses(ipulse+1) = timing(iTrial).motionon + (plength*ipulse);
        end
        
        timing(iTrial).plxstart = syncClock.plexonTrialStartTime(iTrial);
    end
    
    % some final cleanup
    %%%%%%%%%%%%%%%%%%%%%
    trialCnt = trialCnt';
    trialId = trialId';
    trialnumber = trialnumber';
    timing = timing';
    
    if ~isequal(xTrials, nTrials)
        basephase  = basephase(1:xTrials, :, :);
        maskphase  = maskphase(1:xTrials, :, :);
        correct    = correct(1:xTrials, :);
        dirprob    = dirprob(1:xTrials, :);
        goodtrial  = goodtrial(1:xTrials, :);
        pulses     = pulses(1:xTrials, :, :);
        sacTime    = sacTime(1:xTrials);
        targ1XY    = targ1XY(1:xTrials, :);
        targ2XY    = targ2XY(1:xTrials, :);
        targchosen = targchosen(1:xTrials, :);
        targcorrect = targcorrect(1:xTrials, :);
        trialCnt   = trialCnt(1:xTrials, :);
        trialId    = trialId(1:xTrials, :);
        trialnumber = trialnumber(1:xTrials, :);
    end
    
    %repackage
    stim.basephase = basephase;
    stim.condition = condition;
    stim.correct   = correct;
    stim.dirprob        = dirprob;
    stim.exname         = exname;
    stim.eyepos         = eyepos;
    stim.frate          = frate;
    stim.frozenTrialIds = frozenTrialIds;
    stim.gaborXY        = gaborXY;
    stim.goodtrial      = goodtrial;
    stim.luminanceMinMuMax = luminanceMinMuMax;
    stim.maskphase = maskphase;
    stim.nFrames = nFrames;
    stim.nGabors = nGabors;
    stim.nPulses  = nPulses;
    stim.nTrials = nTrials;
    stim.ppd   = ppd;
    stim.pulses = pulses;
    stim.sf = sf;
    stim.targ1XY = targ1XY;
    stim.targ2XY = targ2XY;
    stim.targchosen = targchosen;
    stim.targcorrect = targcorrect;
    stim.tf = tf;
    stim.theta = theta;
    stim.timing = timing;
    stim.trialCnt = trialCnt;
    stim.trialId = trialId;
    stim.trialnumber = trialnumber;
    
    %clear ans iTrial jTrial ipulse iPulse plength tmpcohix trialMode trIx dirDataPDS pdslist PDS trialMatrix iSession syncClock screenLatency
    %save(['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo/stim/' exname '_stim'])
    save(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo/stim/' exname '_stim'], '-struct', 'stim', '-v7.3')    
    disp(['saved exp ' exname '_stim']);
    clearvars -except dirDataPDS pdslist;
end % iSession