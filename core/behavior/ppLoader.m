function [ppAll] = ppLoader(subj, eList)

% Quick wrapper function to load a desired set of pp1 files and place them
% all in a struct.
% INPUTS: subj  - three letter subject code string
%         eList - list of experiments to load (e.g. 1:7)
%         %%filetype - pp, gg, or ss
% OUTPUT: struct containing all your desired pp1 files.

[~, host] = system('hostname');

if strcmp(host(1), 'd') %dhcp-129-116-178-237.cps.utexas.edu
    %desktop
    dirData = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/pp/'; 
else
    %laptop
    dirData = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/pp/';
end

% dirBase     = '/Users/Aaron/Dropbox/_twag/';
% dirData     = fullfile(dirBase, 'data/');

nExp  = length(eList);
count = 1;

% switch filetype
%     case 'pp'

% load desired pp1files
%for ii=eList(1):eList(nExp)
for ii=1:nExp
    if eList(ii) < 10
        ppAll{count, 1} = load([dirData subj '-e0' num2str(eList(ii)) '-pp.mat']);
        %ppAll{count, 1} = load([dirData 'pp1-' subj '-e0' num2str(ii) '.mat']);
    else
        ppAll{count, 1} = load([dirData subj '-e' num2str(eList(ii)) '-pp.mat']);
        %ppAll{count, 1} = load([dirData 'pp1-' subj '-e' num2str(ii) '.mat']);
    end
    disp('Loading data set...');
    count = count+1;
end
disp('Done loading data sets.');

%     case 'gg'
%          for ii=eList(1):eList(nExp)
%             if ii < 10
%                 ppAll{count, 1} = load([dirData '/final/gg1-' subj '-e0' num2str(ii) '.mat']);
%             else
%                 ppAll{count, 1} = load([dirData 'pp1-' subj '-e' num2str(ii) '.mat']);
%             end
%             disp('Loading data set...');
%             count = count+1;
%         end
%         disp('Done loading data sets.');