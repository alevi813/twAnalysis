
% function sacTime = getTwagSaccadeTimes(PDS, screenLatency)

%test session
%load('/Volumes/LKCLAB/EPHYS/DATA/Projects/twag/nancy/20160726/nTwag20160726_mt_l5p4_exp1_final/nancy20160726twag.twag1540.PDS', '-mat')

load('/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy/nancy20160804twag.twag1441.PDS', '-mat'); % laptop file
%% stimulus timings

%if nargin < 2
screenLatency = 0.0524634; % thats 52.4624ms
%end

sacTime    = nan(length(PDS.data), 1);
sacBinDiff = nan(length(PDS.data), 1);

for iTrial = 1:length(PDS.data)
    
    if PDS.data{iTrial}.pldaps.goodtrial ~= 1
        % do nothing
        fpoffTime(iTrial) = NaN;
        sacTime(iTrial)   = NaN;
    else
        
        fpoffTime(iTrial) = PDS.data{iTrial}.stimulus.statesStartTime(6)+ screenLatency;
        
        %% eyelink data & timings
        
        sRate    = PDS.initialParametersMerged.eyelink.srate;
        eyetime = (1:PDS.data{iTrial}.eyelink.sampleNum) / sRate;
        
        %%---- eyelink traces
        eyex = PDS.data{iTrial}.eyelink.samples(4,:);
        eyey = PDS.data{iTrial}.eyelink.samples(5,:);
        
        %%---- round to 3pts and make sure divisible by 2... workaround
        fpRound     = round(fpoffTime(iTrial), 3, 'decimals');
        
        if mod(fpRound, 1/sRate) > 0
            fpRound = fpRound - .001;
        end
        
        % eyelink sample/bin nearest to go signal
        try
            fpoffSample = find(eyetime == fpRound);
        catch ME
        end
        
        if isempty(fpoffSample)
            [~, fpoffSample] = min(abs(eyetime-fpRound));
        end
        
        assert(~isempty(fpoffSample), 'No fpoffsample')
        
        % trial time from go signal
        fromGoTime = eyetime(fpoffSample:end);
        
        % eye trace from go signal, diffed
        fromGoEyex = eyex(fpoffSample:end);
        fromGoEyex = diff(fromGoEyex);
        
        fromGoEyey = eyey(fpoffSample:end);
        fromGoEyey = diff(fromGoEyey);
        
        %         trace = [fromGoEyex; fromGoEyey];
        %         result = saccadeDetector(fromGoTime, trace);
        
        %% ---- find saccade

        %%---- x pos
        [~, maxx_ix] = max(abs(fromGoEyex));
        if maxx_ix < 11
            [~, maxx_ix] = max(abs(fromGoEyex(11:end)));
            
            sacBinx = maxx_ix+10;
        else
            
            sacBinx = maxx_ix;
        end        
        
        %%---- y pos
        [~, maxy_ix] = max(abs(fromGoEyey));
        if maxy_ix < 11
            [~, maxy_ix] = max(abs(fromGoEyey(11:end)));
            
            sacBiny = maxy_ix+10;
        else
            
            sacBiny = maxy_ix;
        end
            
%         sacThresh = 3; % how to optimize this value, hmmmm
%         
%         %         sacBinx = find(fromGoEyex > sacThresh);
%         %         sacBiny = find(fromGoEyey > sacThresh);
%         sacBinx = find(abs(fromGoEyex) > sacThresh);
%         sacBiny = find(abs(fromGoEyey) > sacThresh);
%         
%         % get consecutive threshold crossings (akin to large velocities for
%         % a sustained duration)
%         threshCrossX = diff(find(fromGoEyex > sacThresh));
%         threshCrossY = diff(find(fromGoEyey > sacThresh));
%         
%         sacBinx = sacBinx(threshCrossX==1);
%         sacBiny = sacBiny(threshCrossY==1);
%         
%         % using 10 as an arbitrary cut off for stuff that's too early?
%         %sacBinx = sacBinx(sacBinx > 10);
%         sacBinx = sacBinx(1);
%         %sacBiny = sacBiny(sacBiny > 10);
%         sacBiny = sacBiny(1);
%         
%         if isempty(sacBinx)
%             sacBinx = NaN;
%         else
%             sacBinx = sacBinx(1);
%         end
%         
%         if isempty(sacBiny)
%             sacBiny = NaN;
%         else
%             sacBiny = sacBiny(1);
%         end                
%         
%         % this could be useful to test if you're not getting real saccades.
%         sacBinDiff(iTrial) = sacBinx - sacBiny;
        
        % use the earlier time as the saccade time
        if sacBinx < sacBiny
            sacTime(iTrial) = fromGoTime(sacBinx);
        else
            sacTime(iTrial) = fromGoTime(sacBiny);
        end
        
    end % goodtrial check
end % trial loop

badTrials = find(abs(sacBinDiff > 10));

if ~isempty(badTrials)
    disp( [PDS.initialParametersMerged.session.file(6:13) ' WARNING: This session may have misidentified saccades on ' num2str(length(badTrials)) ' trials '] );
end
