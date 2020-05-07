function [pp] = twagConvertFromCell(pdsf, dirData, trStart)
%% load and convert everything to doubles/matrices for convenient analysis
% Inputs:  pds full path and file name
% Outputs: a smaller struct version of the pds (pp) with necessary
%          fields for later analysis.
%

if nargin < 3
    trStart = NaN;
    if nargin < 2
        dirData = 0;
    end
end
%% load file (from input, pdsf)
if dirData ~= 0
    load(fullfile(dirData, pdsf), '-mat');
else 
    load(pdsf, '-mat');
end
%%

% get your goodtrials straightened out, then cell fun everything
nPulses = cellfun(@(x) x.theseGabs.nPulses, PDS.data, 'UniformOutput', false);

% do a quick silly check for missing fields. These will happen due to
% pauses
for ii = 1:length(PDS.data)
    if ~isfield(PDS.data{ii}, 'targRchosen')
        PDS.data{ii}.targRchosen = nan;
    end
    
    if ~isfield(PDS.data{ii}, 'correct')
        PDS.data{ii}.correct = nan;
    end
    
    if ~isfield(PDS.data{ii}.theseGabs, 'targCorrect')
        PDS.data{ii}.theseGabs.targCorrect = nan;
    end
end

goodtrial = cellfun(@(x) x.pldaps.goodtrial, PDS.data, 'UniformOutput', false);
goodtrial = cell2mat(goodtrial);
goodtrial(isnan(goodtrial)) = 0;

cohmat       = cellfun(@(x) x.theseGabs.cohmat, PDS.data, 'UniformOutput', false);
sumCoh       = cellfun(@(x) x.theseGabs.sumCoh, PDS.data, 'UniformOutput', false);
stimDistNum  = cellfun(@(x) x.theseGabs.stimDistNum, PDS.data, 'UniformOutput', false);
stimDistName = cellfun(@(x) x.theseGabs.stimDistName, PDS.data, 'UniformOutput', false);
targCorrect  = cellfun(@(x) x.theseGabs.targCorrect, PDS.data, 'UniformOutput', false);

correct      = cellfun(@(x) x.correct, PDS.data, 'UniformOutput', false);
targRchosen  = cellfun(@(x) x.targRchosen, PDS.data, 'UniformOutput', false);

% cell2mat everything. Now you should have some workable variables
pp.nPulses       = cell2mat(nPulses(goodtrial==1));
pp.cohmat        = cell2mat(cohmat(goodtrial==1));
       pp.cohmat = pp.cohmat';
pp.sumCoh        = cell2mat(sumCoh(goodtrial==1));
pp.stimDistNum   = cell2mat(stimDistNum(goodtrial==1));
% pp.stimDistName  = cell2mat(stimDistName(goodtrial==1));
pp.targCorrect   = cell2mat(targCorrect(goodtrial==1));
pp.correct       = cell2mat(correct(goodtrial==1));
pp.targRchosen   = cell2mat(targRchosen(goodtrial==1));

pp.stimDistName  = stimDistName(goodtrial==1);
pp.goodtrial     = goodtrial;

pp.initialParametersMerged =  PDS.initialParametersMerged;

% save frozen trial index, if it exists
if isfield(PDS.data{1}.stimulus, 'frozenTrials')
    frozenTrials    = cellfun(@(x) x.stimulus.frozenTrials, PDS.data, 'UniformOutput', false);
    pp.frozenTrials = cell2mat(frozenTrials(goodtrial==1));
end

if isfield(PDS.data{1}, 'contrast')
    contrast   = cellfun(@(x) x.contrast, PDS.data, 'UniformOutput', false);
    pp.contrast = cell2mat(contrast(goodtrial==1));
    pp.contrast = pp.contrast*2;
end

%%

if ~isnan(trStart)
pp.nPulses       = pp.nPulses(trStart:end);
pp.cohmat        = pp.cohmat(trStart:end, :);
pp.sumCoh        = pp.sumCoh(trStart:end);
pp.stimDistNum   = pp.stimDistNum(trStart:end);
pp.targCorrect   = pp.targCorrect(trStart:end);
pp.correct       = pp.correct(trStart:end);
pp.targRchosen   = pp.targRchosen(trStart:end);

pp.stimDistName  = pp.stimDistName(trStart:end);
end

