function popcoh = conditionAvgPTA(S, condition, spNum)

%%

field='ptaRaw';
nNeurons = length(S);

% set up color
switch condition
    case 'flat'
        clr = [19/255 40/255 230/255];
    case 'late'
        clr = [255/255 234/255 0/255];
    case 'early'
        clr = [214/255 4/255 0/255];
end
%%% load data from S
% normalize by max
cohpData = reshape(cell2mat(arrayfun(@(x) x.model(1).(field)/max(x.model(1).(field)(:)), S, 'UniformOutput', false)), [size(S(1).model(1).(field)) nNeurons]);

% don't normalize by max
%cohpData = reshape(cell2mat(arrayfun(@(x) x.model(1).(field), S, 'UniformOutput', false)), [size(S(1).model(1).(field)) nNeurons]);

popcoh   = nanmean(cohpData,3);

% set up time axis
timex=S(1).model(1).ptaTime;
xdim = [0 2000];

timex = timex(1:50, :);
popcoh = popcoh(1:50, :);
xdim = [0 1400];

% plot-
subplot(3,3, spNum)
for kP=1:size(popcoh,2)
    if size(timex,2)>1
        txx=timex(:,kP);
    else
        txx=timex;
    end
        
    plot(txx, popcoh(:,kP), 'Color', clr, 'LineWidth', 2); hold on
end

xlim(xdim)
axis square
set(gcf, 'color', 'white');

