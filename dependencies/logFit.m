function [xfit, fitvals] = logFit(xvals, yvals)

xfit = linspace(xvals(1), xvals(end));

logfun = @(x, xvals) (x(1)+x(2)*log(xvals)); x0 = [0, .5];

xhat = lsqcurvefit(logfun,x0,xvals,yvals);

fitvals = logfun(xhat,xfit);
