function [xvals, respFit] = dsFit(thetas, responses)

opt=optimset;
startvals = [5  5       pi       1     1];
ub =        [99  99   2*pi      10    35]; 
lb =        [1   1       0       0     0]; 


[x,rn,res,ef,out,lam,jac] = lsqcurvefit('fitgauss',startvals,thetas.*(pi/180),responses',lb,ub,opt);

% if cc==26
%     x = [3.6422    2   1.6    0.2334    2];
% end

xvals=(0:359)';
respFit=fitgauss(x,(0:359).*(pi/180))';

% plot(xvals(1:330),respFit(1:330),'color',[.5 .5 .5],'linewidth',2)