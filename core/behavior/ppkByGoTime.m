
nBins = 6;

goTime = late.goTime-late.motionOn;
h=histogram(goTime, nBins);

figure
for kBin = 1:nBins
    % get bin indicies
    bin(kBin).ix = goTime> h.BinEdges(kBin) & goTime <h.BinEdges(kBin+1); 
    
    % calculate ppk
    ppk = ppkTools(late.pulses(bin(kBin).ix, :), late.cho(bin(kBin).ix, :));
    
    % plot
    subplot(1,nBins, kBin)
    plot(ppk);
    title(num2str(mean([h.BinEdges(kBin) h.BinEdges(kBin+1)])))
end

supertitle('time to "go"', 12)