function spikeData = readRez(rez)
%function spikeData = readRez(rez)

sampleRate = rez.ops.fs;
spikeSamples = rez.st3(:,1);
spikeTimes = spikeSamples ./ sampleRate; % convert samples to times!
clusters = rez.st3(:,5); %+ 1; % add 1? since they start at zero (b/c of python indexing) = whatever let's leave it off
clusterIDs = unique(clusters);
spikeTemplates = rez.st3(:,2);
amplitudes = rez.st3(:,3);

%% get channelIDs corresponding to clusters based on peak amplitude
% % peakChannels = zeros(size(clusters));
% % for c = 1:length(clusterIDs)
% %     clust     = clusterIDs(c);
% %     I         = clusters == clust;
% %     templates =  unique(spikeTemplates(I));
% %     
% %     t = squeeze(range((rez.dWU(:,:,templates)),1));
% %     m = max(max(t));
% %     
% %     try
% %         if all(m)
% %             if any(size(t) == 1)
% %                 chidx = find(t == m);
% %             else
% %                 [chidx, ~] = find(t == m);
% %             end
% %             peakChannels(I) = chidx; % *** some sessions this breaks - not sure why, that's why I'm adding a try... catch
% %         else % if m is zero (for debugging)
% %             peakChannels(I) = NaN;
% %         end
% %     catch
% %         fprintf('\n *** Peak channel indexing issue, be aware. *** \n')
% %         peakChannels = nan(size(clusters));
% %     end
% % end

% this doesn't even put it back by each unit (cluster ID) anymore... look
% at this if you go back to this

%% new solution that I was using when getting the stuff directly from the the py files, not sure if this will work here jsut copying and pasting because the indexing might be off, or variables might now be eactly the same, we'll see
temps = permute(rez.Wraw,[3 2 1]); % looks like it was there, but in a different order, so have to rearrange

% these were from the function to read in the kilosort py files
% clusters = spikeStruct.clu;
% clusterIDs = spikeStruct.cids;
% spikeTemplates = spikeStruct.spikeTemplates;  % note: zero-indexed !!!!!!!!!!!!!!!!!!!
% amplitudes = spikeStruct.tempScalingAmps;
%temps = spikeStruct.temps; % [nTemplates, nTimePoints, nTempChannels]

peakChannels = nan(1,length(clusterIDs));
for c = 1:length(clusterIDs)
    clust     = clusterIDs(c);
    I         = clusters == clust;
    templates =  unique(spikeTemplates(I));%+1; % +1 because it was zero indexed!!!! (that was only if it was coming directly from the py file I believe)
    
    % temps = [nTemplates, nTimePoints, nTempChannels]
    t = squeeze(range((temps(templates,:,:)),1)); % what dimension to range over???????
    
    maxWFperCh = max(t);
    maxWF = max(max(t));
    if ~all(maxWFperCh==0)
        peakCh = find(maxWFperCh==maxWF);
    else
        peakCh = NaN;
    end
    
    peakChannels(c) = peakCh;
end

% so ... can't tell if this is any better or worse... have to look into it

%%
spikeData.sampleRate = sampleRate;
spikeData.spikeTimes = spikeTimes;
spikeData.clusters = clusters;
spikeData.clusterIDs = clusterIDs;
spikeData.peakChannel = peakChannels;
spikeData.spikeSamples = spikeSamples;
spikeData.spikeTemplates = spikeTemplates;
spikeData.amplitudes = amplitudes;

% NOTE some code was you using peakChannel and others peakChannels OMFG
% erghhhh, changing everything to peakChannels