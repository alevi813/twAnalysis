function [cohBins] = binCoherences(X, nbins)
% bin stimulus strength with different spacings
% [coh, abscoh, binCenters, binEdges, absBinCenters, absBinEdges] = binCoherences(X, nbins, binopts(not used yet))

if size(X,2) == 1
    coh = X;
else
    coh     = sum(X,2); % signed
end
abscoh  = abs(coh); % absolute

% -------------------------------------------------------------------------
% bin coherence levels
[~,binCenters] = hist(coh, nbins);
%             binCenters = [-.6 -.5 -.3 -.2 -.1 .1 .2 .3 .4 .5];
%             nbins = numel(binCenters);
% add bins around zero
% zbin       = find(binCenters.^2 < 1e-3); % bin centered at zero
% binCenters = [binCenters(1:zbin-1) binCenters(zbin-1:zbin+1)/2 binCenters(zbin+1:end)];
binEdges   = [binCenters(1:2) - diff(binCenters(1:2))/2 binCenters(2:end) + diff(binCenters)/2];
% absolute coherence bins
absBinEdges   = [0 binEdges(binEdges >=0)];
absBinCenters = absBinEdges(1:end-1) + diff(absBinEdges)/2;

% reshape into column vectors
coh = coh(:);
abscoh = abscoh(:);
binCenters = binCenters(:);
binEdges = binEdges(:);
absBinCenters = absBinCenters(:);
absBinEdges = absBinEdges(:);

id = zeros(size(coh));
for k = 1:numel(binCenters)
    ndx = coh > binEdges(k) & coh < binEdges(k+1);
    id(ndx) = k;
end

cohBins.coh = coh;
cohBins.abscoh = abscoh;
cohBins.binCenters = binCenters;
cohBins.binEdges = binEdges;
cohBins.absBinCenters = absBinCenters;
cohBins.absBinEdges = absBinEdges;
cohBins.id = id;

end