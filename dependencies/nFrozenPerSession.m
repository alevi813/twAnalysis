leoDir    = ('/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo/stim/');
stimFiles = dir(leoDir);

stimFiles = stimFiles(~ismember({stimFiles.name},{'.','..','.DS_Store'})); % CLEAN UP
nLeoFiles = length(stimFiles);

for iS = 1:nLeoFiles
   
    stim  = load([leoDir stimFiles(iS).name]);
    froIx = stim.trialId==stim.frozenTrialIds;
    froIx = froIx(stim.goodtrial);
    
    nFro(iS) = sum(froIx);
end

%%

nancyDir    = ('/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/nancy/stim/');
stimFiles = dir(nancyDir);

stimFiles = stimFiles(~ismember({stimFiles.name},{'.','..','.DS_Store'})); % CLEAN UP
nNancyFiles = length(stimFiles);

for iS = 1:nNancyFiles
   
    stim  = load([nancyDir stimFiles(iS).name]);
    froIx = stim.trialId==stim.frozenTrialIds;
    froIx = froIx(stim.goodtrial);
    
    nFro(iS+nLeoFiles) = sum(froIx);
end