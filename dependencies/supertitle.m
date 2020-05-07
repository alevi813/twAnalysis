function [] = supertitle(titleIn, fontSize)

% the function uses the annotation function to place a title at the top
% part of a figure.
% 
% 2013 lnk wrote it
% 20140522 lnk added a fontSize input variable

if nargin < 2
    fontSize = 6;
end


annotation('textbox', [.02 .98 .98 .02], 'string', titleIn, 'linestyle', 'none', 'fontsize', fontSize)
        