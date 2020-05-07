function [b, dev, s]=computeKernel(cohmat, distIdx, targRchosen, normalize)

%% compute kernel

cohmat   = cohmat(distIdx==3,:);
cohmat   = cohmat/std(cohmat(:));
choice   = targRchosen(distIdx==3);

[b,dev,s]=glmfit(cohmat,choice','binomial');
%[b,dev,s]=glmfit(cohmat,choice,'binomial');

if normalize
    b  = b/norm(b);
    s.se = s.se/norm(b);
end

% toggle this on/off if you don't need/do need the bias term
b = b(2:end);
s.se = s.se(2:end);


