function S = loadNeuronsByCondition(condition, monkey)

%% load all neurons

if nargin < 2
    monkey = 'both';
end

comp = getComp;

% path info
if strcmp(comp, 'laptop')
    dataPath{1} = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/nancy';
    dataPath{2} = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo';
else
    dataPath{1} = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy';
    dataPath{2} = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo';
end

% load nancy neurons
if ~strcmp(monkey, 'leo')
    fitDir    = 'fits_cho1'; % main_fits / fits_notargs / fits_pulsecovars / fits_cho1
    S_nancy = loadUpFits(dataPath{1}, 'MT', condition, fitDir); % rate predictions under choice 1 kernel (right)
end

% load leo neurons
if ~strcmp(monkey, 'nancy')
    fitDir    = 'fits_leo'; % main_fits / fits_notargs / fits_pulsecovars / fits_cho1
    S_leo = loadUpFits(dataPath{2}, 'MT', condition, fitDir); % rate predictions under choice 1 kernel (right)
end


% combine
if strcmp(monkey, 'both')
    S = [S_nancy S_leo];
else
    if exist('S_nancy', 'var')
        S = S_nancy;
    end
    
    if exist('S_leo', 'var')
        S = S_leo;
    end
end