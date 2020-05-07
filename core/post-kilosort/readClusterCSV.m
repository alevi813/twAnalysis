function [clusterID, clusterGroup] = readClusterCSV(exDir)

%% Read in manually sorted cluster info
filename = [exDir filesep 'Kilosort' filesep 'cluster_groups.csv'];
delimiter = '\t';
startRow = 2;
formatSpec = '%f%s%[^\n\r]'; %   column1: double (%f), column2: text (%s)
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
clusterID = dataArray{:, 1};
clusterGroup = dataArray{:, 2};