function cp = twagSessionCP(cpPlot, cp, goodChannelIx, bc)
% calculate grand mean CP for the whole session and then pulse CP for whole session
% for all revco and frozen trials

%% all revco
cp.session.allRevco.m = cell2mat(cp.channel.allRevco.m(goodChannelIx));
cp.session.allRevco.m = reshape(cp.session.allRevco.m, [1700, sum(goodChannelIx)]);
cp.session.allRevco.m = mean(cp.session.allRevco.m, 2);

cp.session.allRevco.s = cell2mat(cp.channel.allRevco.s(goodChannelIx));
cp.session.allRevco.s = reshape(cp.session.allRevco.s, [1700, sum(goodChannelIx)]);
cp.session.allRevco.s = mean(cp.session.allRevco.s, 2);

if cpPlot == 1
    figure; hold on
    plot(bc, cp.session.allRevco.m,'k',...
        bc, cp.session.allRevco.m+cp.session.allRevco.s, 'k--',...
        bc, cp.session.allRevco.m-cp.session.allRevco.s, 'k--');
    xlim([0 1.5]);
    
    xlabel('time')
    ylabel('CP')
    title('Session average All Revco')
end

% % pulse CP all revCo
% for ii = 1:7
%     sessionPulse_mCP_allRevco(ii, :) = cp.session.allRevco.m(pStartBins(ii):(pStartBins(ii)+pLength));
%     %sessionPulse_sCP(ii, :) = session_sCP(pStartBins(ii):(pStartBins(ii)+pLength));
% end
% 
% cp.session.allRevco.pulse.m    = mean(sessionPulse_mCP_allRevco,2);
% cp.session.allRevco.pulse.std  = std(sessionPulse_mCP_allRevco');
% cp.session.allRevco.pulse.se   = sessionPulseCPstd_allRevco/sqrt(length(sessionPulse_mCP_allRevco));

%% frozen only
cp.session.frozen.m = cell2mat(cp.channel.frozen.m(goodChannelIx));
cp.session.frozen.m = reshape(cp.session.frozen.m, [1700, sum(goodChannelIx)]);
cp.session.frozen.m = mean(cp.session.frozen.m, 2);

cp.session.frozen.s = cell2mat(cp.channel.frozen.s(goodChannelIx));
cp.session.frozen.s = reshape(cp.session.frozen.s, [1700, sum(goodChannelIx)]);
cp.session.frozen.s = mean(cp.session.frozen.s, 2);

if cpPlot == 1
    figure; hold on
    plot(bc, cp.session.frozen.m,'k',...
        bc, cp.session.frozen.m+cp.session.frozen.s, 'k--',...
        bc, cp.session.frozen.m-cp.session.frozen.s, 'k--');
    xlim([0 1.5]);
    
    xlabel('time')
    ylabel('CP')
    title('Session average Frozen only')
end

% % pulse CP frozen
% for ii = 1:7
%     sessionPulse_mCP_frozen(ii, :) = cp.session.frozen.m(pStartBins(ii):(pStartBins(ii)+pLength));
%     %sessionPulse_sCP(ii, :) = session_sCP(pStartBins(ii):(pStartBins(ii)+pLength));
% end
% 
% cp.session.frozen.pulse.m    = mean(sessionPulse_mCP_frozen,2);
% cp.session.frozen.pulse.std  = std(sessionPulse_mCP_frozen');
% cp.session.frozen.pulse.se   = sessionPulseCPstd_frozen/sqrt(length(sessionPulse_mCP_frozen));
