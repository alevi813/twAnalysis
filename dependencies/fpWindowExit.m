function exitSample = fpWindowExit(trace, winDims, sRate)
% give the trace in a 2xnSample matrix. Row 1 is xtrace, row 2 is ytrace
% winDims are pixel dimensions of the fixation window. 1&3 are x. 2&4 are y


x_outOfWindow = trace(1,:) < winDims(1) | trace(1,:) > winDims(3);
y_outOfWindow = trace(2,:) < winDims(2) | trace(2,:) > winDims(4);

x_crossIx = find(x_outOfWindow);
y_crossIx = find(y_outOfWindow);

if isempty(x_crossIx)
    x_crossIx = Inf;
end
if isempty(y_crossIx)
    y_crossIx = Inf;
end

if x_crossIx(1) < y_crossIx(1)
    exitSample = x_crossIx(1)-(50*(sRate*1e-3));
    if exitSample < 1
        exitSample = x_crossIx(1);
    end    
else
    exitSample = y_crossIx(1)-(50*(sRate*1e-3));
    if exitSample < 1
        exitSample = y_crossIx(1);
    end
end

    
