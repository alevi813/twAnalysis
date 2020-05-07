% ptaPostKilo
subject = 'leo';
sessionList = 1:10;
dataDir = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData';

for isession = 1:length(sessionList)
    
   session = load(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/kiloSortAnalysis/processedData' filesep subject 'Session' num2str(sessionList(isession)) '.mat']);
   
   [sessionPTA, sessionMaxPTA] = unitPTA_postKilo(session, 1, sessionList(isession)); 
    
end