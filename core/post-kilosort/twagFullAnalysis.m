%% analyses using grand avg of all neurons
%  cleaning things up and centralizing analyses here for final figures

comp = getComp;

psth      = false;
cp        = true;   shuff = false;
pta       = false;
noisecorr = false;
nmf       = false;
dpHist    = false;
%% load all neurons, loop through conditions plotting psth, cp, pta


condition = {'flat', 'late', 'early'}; % flat/late/early --- use 'good' for all sessions
%condition = {'late'}; % flat/late/early --- use 'good' for all sessions

for iCond = 1:length(condition)
    
    switch condition{iCond}
        case 'flat'
            clr = [19/255 40/255 230/255];
        case 'late'
            clr = [255/255 234/255 0/255];
        case 'early'
            clr = [214/255 4/255 0/255];
        case 'good'
            clr = [.45 .45 .45];
    end
    
    % path info
    if strcmp(comp, 'laptop')
        dataPath{1} = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/nancy';
        dataPath{2} = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo';
    else
        dataPath{1} = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy';        
        dataPath{2} = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo';
    end
        
    % load nancy neurons
    fitDir    = 'fits_cho1'; % main_fits / fits_notargs / fits_pulsecovars / fits_cho1
    S_nancy = loadUpFits(dataPath{1}, 'MT', condition{iCond}, fitDir); % rate predictions under choice 1 kernel (right)
    
    % load leo neurons
    fitDir    = 'fits_leo'; % main_fits / fits_notargs / fits_pulsecovars / fits_cho1
    S_leo = loadUpFits(dataPath{2}, 'MT', condition{iCond}, fitDir); % rate predictions under choice 1 kernel (right)
    
    % threshold out cells w/ low d'
    dpThresh = false; % doing this later now instead of at the beginning
    if dpThresh
        dp = arrayfun(@(x) x.model(1).dprime, S_nancy);
        dpix = abs(dp) >= 0.1; %0.15
        S_nancy = S_nancy(dpix);
        
        dp = arrayfun(@(x) x.model(1).dprime, S_leo);
        dpix = abs(dp) >= 0.1; %0.15
        S_leo = S_leo(dpix);
    end
    
    % combine
    S = [S_nancy S_leo];
    
    dpAll = arrayfun(@(x) x.model(1).dprime, S);
    dpNancy = arrayfun(@(x) x.model(1).dprime, S_nancy);
    dpLeo = arrayfun(@(x) x.model(1).dprime, S_leo);
    
    if dpHist
        figure(1); 
        spNum = [2 3 1];
        subplot(3,3,spNum(iCond)); hold on
        histogram(abs(dpAll), 'FaceColor', clr)
        title(['med abs(dp) = ' num2str(nanmedian(abs(dpAll)))])
        axis square

        spNum = [5 6 4];
        subplot(3,3,spNum(iCond))
        histogram(abs(dpLeo), 'FaceColor', clr)
        title(['med abs(dp) = ' num2str(nanmedian(abs(dpLeo)))])
        axis square

        spNum = [8 9 7];
        subplot(3,3,spNum(iCond))
        histogram(abs(dpNancy), 6, 'FaceColor', clr)
        title(['med abs(dp) = ' num2str(nanmedian(abs(dpNancy)))])  
        axis square

    end
    
    %% PSTH
    if psth
        % plot psth by stimulus distribution
        spNum = [2 3 1];
        figure(1)
        responseByDist(S, dataPath, spNum(iCond));
        
        spNum = [5 6 4];
        responseByDist(S_leo, dataPath, spNum(iCond));
        
        spNum = [8 9 7];
        responseByDist(S_nancy, dataPath, spNum(iCond));
    end
    
    %% Choice Probability
    if cp
        figure(2)
        spNum = [2 3 1];
        [cpm, cpse, psthTime] = plotCPt(S, dataPath, spNum(iCond));
        %pulseCP(S, dataPath, spNum(iCond));
        
        figure(3)
         subplot(3,3,spNum(iCond))
         cpAll = grandCP_test(S, dataPath, shuff);
        h = histogram(cpAll, 'FaceColor', clr);
        h.BinWidth = 0.05; 
        title(['med cp = ' num2str(nanmedian(cpAll))])  
        axis square
        xlim([0.05 .95])
        
        figure(2)
       spNum = [5 6 4];
       plotCPt(S_leo, dataPath, spNum(iCond));
       % pulseCP(S_leo, dataPath, spNum(iCond));
        
       figure(3)
         subplot(3,3,spNum(iCond))
         cpLeo = grandCP_test(S_leo, dataPath, shuff);
        h=histogram(cpLeo, 'FaceColor', clr);
        h.BinWidth = 0.05; 
        title(['med cp = ' num2str(nanmedian(cpLeo))])  
        axis square
        xlim([0.05 .95])

        figure(2)
       spNum = [8 9 7];
       plotCPt(S_nancy, dataPath, spNum(iCond));
       % pulseCP(S_leo, dataPath, spNum(iCond));
        
       figure(3)
         subplot(3,3,spNum(iCond))
         cpNancy = grandCP_test(S_nancy, dataPath, shuff);
         h=histogram(cpNancy, 'FaceColor', clr);
         h.BinWidth = 0.05; 
         title(['med cp = ' num2str(nanmedian(cpNancy))])  
         axis square
         xlim([0.05 .95])

        if shuff
            % get a 2sem range
            % all neurons --- flat:  0.0088;  late:  0.0079;   early: 0.0086;
            ci_all(iCond)   = 2* (std(cpAll) / sqrt(length(cpAll)) ) ;
            ci_leo(iCond)   = 2* (std(cpLeo)  / sqrt(length(cpLeo)) ) ;
            ci_nancy(iCond) = 2* (std(cpNancy) / sqrt(length(cpNancy)) ) ;

            %cpSig = cpAll > (mean(cpAll)+ ci) | cpAll < (mean(cpAll)- ci);
        end

    end
    
    %% Pulse-triggered average
    if pta
        figure(13)
        spNum = [2 3 1];
        pta_data{iCond} = conditionAvgPTA(S, condition{iCond}, spNum(iCond));
        
        spNum = [5 6 4];
        conditionAvgPTA(S_leo, condition{iCond}, spNum(iCond));
        
        spNum = [8 9 7];
        conditionAvgPTA(S_nancy, condition{iCond}, spNum(iCond));
    end
    
    %% dp vs. cp    
%     figure(4)
%     
%     cpAll = abs(cpAll - .5); %%%% test
%     spNum = [2 3 1];
%     subplot(3,3,spNum(iCond))
%     plot(abs(dpAll), cpAll, 'o', 'Color', clr)
%     r = corrcoef(abs(dpAll), cpAll); %%%% r = 0.037, all units.
%     title([' r = ' num2str(r(1,2))])
%     axis square
%     
%     cpLeo = abs(cpLeo - .5); %%%% test
%     spNum = [5 6 4];
%     subplot(3,3,spNum(iCond))
%     plot(abs(dpLeo), cpLeo, 'o', 'Color', clr)
%     r = corrcoef(abs(dpLeo), cpLeo);
%     title([' r = ' num2str(r(1,2))])
%     axis square
%     
%     cpNancy = abs(cpNancy - .5); %%%% test
%     spNum = [8 9 7];
%     subplot(3,3,spNum(iCond))
%     plot(abs(dpNancy), cpNancy, 'o', 'Color', clr)
%     r = corrcoef(abs(dpNancy), cpNancy);
%     title([' r = ' num2str(r(1,2))])
%     axis square    
    
    %% interneuronal correlations
    if noisecorr
    %figure(5)
    spNum = [2 3 1];
    [rval, tmpR_means] = interNeuronalCorrelations(S, dataPath, spNum(iCond));
    rMeans(spNum(iCond)) = tmpR_means;
%     plot(median(cell2mat(rval.first(:))), 225, 'bv', 'markersize', 15)
%     plot(median(cell2mat(rval.last(:))), 225, 'rv', 'markersize', 15)
    end
    
    %% neurometric functions
    if nmf
        [xAll, allAuc, ~] = nmfByNeuron(S, dataPath);
        
        errorbar(nanmean(xAll), nanmean(allAuc), nanstd(allAuc)/sqrt(size(allAuc, 1)) ,'d', 'color', [.6 .6 .6])
        [xfit, fitvals] = logFit(nanmean(xAll), nanmean(allAuc));
        hold on
        plot(xfit, fitvals, 'color', [.6 .6 .6])

    end
end




