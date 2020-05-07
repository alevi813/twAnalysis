function [sessionPTA, sessionMaxPTA] = unitPTA_postKilo(session, ptaPlot, nsession)

if nargin < 3
    ptaPlot = 0;
end
%figDir = '/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/figures';
figDir = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/figures';

for iUnit = 1:length(session.goodClustIx)
    
    theseSpikes = double(session.spikeTimesPtb(session.clusters == session.clusterID(iUnit)));
    
    if session.goodClustIx(iUnit)
        % be careful of all trial vs rc only confusion w/respect to
        % centerTime and cohmat
        [PTA, maxPTA] = pulseTriggeredAverage(ptaPlot, session.pp.cohmat, session.pp.stimDistNum, theseSpikes, session.centerTimeGT, session.dirPref(iUnit), 0);
        
        if ptaPlot
            set(figure(1), 'PaperSize', [16 12], 'PaperPosition', [0 0 16 12])
            saveas(figure(1), [figDir '/leo/ptaAllTrial/pta_session' num2str(nsession) '_unit' num2str(iUnit)], 'epsc');
            %saveas(figure(1), [figDir '/ptaRConly/pta_session' num2str(nsession) '_unit' num2str(iUnit)], 'epsc');
            %saveas(figure(1), ['/Users/Aaron/Desktop/pta_test/' num2str(nsession) '_unit' num2str(iUnit)], 'epsc');
            disp(['saved pta (FIGS) session ' num2str(nsession)]);
            
            close all
        end

    sessionPTA{iUnit}       = PTA;
    sessionMaxPTA(iUnit, :) = maxPTA;     
    end % if goodclust
   
end % for all units

save(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/PTA/leoSession' num2str(nsession)]);
%save(['/Users/Aaron/Dropbox/twagAnalysis4.1/kiloSortAnalysis/PTA/rcOnly/nancySession' num2str(nsession)]);
%save(['/Users/Aaron/Desktop/pta_test/' num2str(nsession)]);
disp(['saved pta (DATA) session ' num2str(nsession)]);