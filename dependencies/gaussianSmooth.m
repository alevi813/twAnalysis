function [smoothed, halfWidth] = gaussianSmooth(data, gausWidth)

gausWidth = int16(gausWidth);
halfWidth = gausWidth/2;
gaussFilter = gausswin(double(gausWidth));
gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
smoothed = conv(data, gaussFilter);