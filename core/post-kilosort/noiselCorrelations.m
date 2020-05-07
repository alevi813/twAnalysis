function [rval, signVal, pairDiff, pairMean] = noiselCorrelations(S, all, dataPath, doShuffle)

% get all experiment/session names
exnames = arrayfun(@(x) x.exname, S, 'UniformOutput', false);
exnames = unique(exnames);

for iExp = 1:length(exnames)
    
    % load the stim file to get dirprob/froIx
    if strcmp(exnames{iExp}(1), 'l')
        stim = load([dataPath{2} filesep 'stim' filesep exnames{iExp} '_stim.mat']);
    else
        stim = load([dataPath{1} filesep 'stim' filesep exnames{iExp} '_stim.mat']);        
    end
    
    froIx = stim.trialId == stim.frozenTrialIds;
    
    if isfield(stim, 'validTrials')
        froIx = froIx(stim.validTrials, :);
    else
        froIx = froIx(stim.goodtrial, :);
    end
    
    % get index present during specific sessoin
    for iNeuron = 1:length(S)
        sessionCells(iNeuron) = strcmp(exnames{iExp}, S(iNeuron).exname);
    end
    
    if sum(sessionCells) > 1
        % get all combinations for units present during the session
        theseCells = S(sessionCells);
        dp = all.dp(sessionCells);
        cp = all.cp(sessionCells);

        allCombs = combnk(1:length(theseCells), 2);
        
        % calculate r_sc for each combination
        for icomb = 1:size(allCombs,1)
            % check if this pair (comb) has the same sign of d' and cp
            signVal.dp{iExp, icomb} = double( isequal( sign(dp(allCombs(icomb, 1))), sign(dp(allCombs(icomb, 2))) ) );
            signVal.cp{iExp, icomb} = double( isequal( cp(allCombs(icomb, 1))>.5, cp(allCombs(icomb, 2))>.5 ) );
            
            if signVal.dp{iExp, icomb} == 1
                % save difference in cp and dp between pairs
                pairDiff.dp{iExp, icomb} = abs(dp(allCombs(icomb, 1))) - abs(dp(allCombs(icomb, 2)));
                % save mean cp and dp between pairs
                pairMean.dp{iExp, icomb} = abs( mean([dp(allCombs(icomb, 1)) dp(allCombs(icomb, 2))]) );
                
                signVal.both{iExp, icomb} = double( isequal( cp(allCombs(icomb, 1))>.5, cp(allCombs(icomb, 2))>.5 ) );
            else
                pairDiff.dp{iExp, icomb} = NaN;
                pairMean.dp{iExp, icomb} = NaN;
                
                signVal.both{iExp, icomb} = NaN;
            end
            
            if signVal.cp{iExp, icomb} == 1
                signVal.cpPos{iExp, icomb} = double( cp(allCombs(icomb, 1))>.5 & cp(allCombs(icomb, 2))>.5 );
                signVal.cpNeg{iExp, icomb} = double( cp(allCombs(icomb, 1))<.5 & cp(allCombs(icomb, 2))<.5 );
                
                tmpcp1 = abs(cp(allCombs(icomb, 1)) - .5);
                tmpcp2 = abs(cp(allCombs(icomb, 2)) - .5);

                % save difference in cp and dp between pairs
                pairDiff.cp{iExp, icomb} = tmpcp1 - tmpcp2;
                % save mean cp and dp between pairs
                pairMean.cp{iExp, icomb} = mean([tmpcp1 tmpcp2]);
            else
                signVal.cpPos{iExp, icomb} = NaN;
                signVal.cpNeg{iExp, icomb} = NaN;
                
                pairDiff.cp{iExp, icomb} = NaN;
                pairMean.cp{iExp, icomb} = NaN;
            end
            
            % let's reassign some stuff to temp vars to make things
            % look a bit cleaner
%             cell1.all    = sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, :), 2);
%             cell1.first  = sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, 100:160), 2);
%             cell1.last   = sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, 160:220), 2);
%             cell1.pre    = sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, 39:99), 2);
%             cell1.post   = sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, 221:280), 2);
% 
%             cell2.all    = sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, :), 2);
%             cell2.first  = sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, 100:160), 2);
%             cell2.last   = sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, 160:220), 2);
%             cell2.pre    = sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, 39:99), 2);
%             cell2.post   = sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, 221:280), 2);
            cell1.all    = zscore(sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, :), 2));
            cell1.first  = zscore(sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, 100:160), 2));
            cell1.last   = zscore(sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, 160:220), 2));
            cell1.pre    = zscore(sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, 39:99), 2));
            cell1.post   = zscore(sum(theseCells(allCombs(icomb, 1)).model(1).trialSpikes(froIx, 221:280), 2));

            cell2.all    = zscore(sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, :), 2));
            cell2.first  = zscore(sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, 100:160), 2));
            cell2.last   = zscore(sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, 160:220), 2));
            cell2.pre    = zscore(sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, 39:99), 2));
            cell2.post   = zscore(sum(theseCells(allCombs(icomb, 2)).model(1).trialSpikes(froIx, 221:280), 2));


            % remove nan's if you've got em
            nanIx1 = isnan(cell1.all); nanIx1 = nanIx1(:,1);
            nanIx2 = isnan(cell2.all); nanIx2 = nanIx2(:,1);
            nanIx = nanIx1+nanIx2;
            nanIx = nanIx > 0;
            
            cell1.all(nanIx, :)   = [];
            cell1.first(nanIx, :) = [];
            cell1.last(nanIx, :)  = [];
            cell1.pre(nanIx, :)   = [];
            cell1.post(nanIx, :)  = [];
            
            cell2.all(nanIx, :)   = [];
            cell2.first(nanIx, :) = [];
            cell2.last(nanIx, :)  = [];
            cell2.pre(nanIx, :)   = [];
            cell2.post(nanIx, :)  = [];
            
            if doShuffle & ~isempty(cell1.all)
                cell1.all   = Shuffle(cell1.all);
                cell1.first = Shuffle(cell1.first);
                cell1.last  = Shuffle(cell1.last);
                
                cell2.all   = Shuffle(cell2.all);
                cell2.first = Shuffle(cell2.first);
                cell2.last  = Shuffle(cell2.last);
            end                

            % calculate r_sc --- full trial
            tempR = corrcoef(cell1.all, cell2.all);
            rval.total{iExp, icomb}  = tempR(1,2);
            
            % calculate r_sc --- first half and second half separately
            tempRfirst = corrcoef(cell1.first, cell2.first);
            tempRlast  = corrcoef(cell1.last, cell2.last);
            rval.first{iExp, icomb}  = tempRfirst(1,2);
            rval.last{iExp, icomb}   = tempRlast(1,2);
            
            % calculate r_sc --- pre and post stimulus
            tempRpre   = corrcoef(cell1.pre, cell2.pre);
            tempRpost  = corrcoef(cell1.post, cell2.post);
            rval.pre{iExp, icomb}    = tempRpre(1,2);
            rval.post{iExp, icomb}   = tempRpost(1,2);            
            
            clear cell1 cell2 
        end % all combinations of neurons
    end % if more than one cell
end % iExp
