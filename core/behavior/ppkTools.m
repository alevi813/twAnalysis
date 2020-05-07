classdef ppkTools < handle
    % Psychophysical Kernel toolbox
    % p = PPKTOOLS(X, Y) - initialize PPKTOOLS object
    %   X is a vector (or matrix) of stimulus strengths, Y is a binary
    %   vector of responses.
    %
    % p = PPKTOOLS(X,Y,'PARAM1', val1, 'PARAM2', val2, ...)
    % 'meanStimulus' - logical to use the mean of the stimulus matrix
    %               instead of each column for the design matrix.
    %
    % p.plot - plots the data and the fit
    % p.plot('PARAM1', val1, ...)
    %   'normWeights'  - plot the normalized weights+se (default=false)
    %   'Color'        - Color of the data and fit
    %   'Marker'       - STRING to set the plot style (default = 'o')
    %
    % plotDiagnostics - gives user a graphic view that's helpful in
    %   troubleshooting. plots data + fit + design matrix + responses.
    %
    % p.xValidate('PARAM1', val1, ...)
    % TBD TBD TBD
    %
    
    properties
        nTrials         % number of trials
        response        % input response values ([nTrials x 1])
        stimulus        % input stimulus values (eg. [nTrials x nPulses])
        designMatrix    % the design matric: zScored stim + bias column
        idx_stimulus    % indices for the stimulus columns in the design matrix
        idx_bias        % indices for the bias column in the design matrix
        costFunction    % cost function used to compute regression weights
        w               % regression weights for design matrix columns
        w_se            % standard error for each weight
        w_norm          % weights, normalized by l2norm
        w_se_norm       % standard error, normalized by l2norm
        w_index         % An index of how early vs late the weighting is
        w_lm            % linear model fit to the weights
        w_com           % center of mass for weights
        w_energy        % sum of sqaured error of residuals from mean of weights
        xvStruct        % struct for the xValidate output
        testZone        % struct for testing output
        hD              % handle for plotted Data
        hF              % handle for plotted Fit
    end
    
    
    
    %%
    methods
        
        %%
        function self = ppkTools(stimulus, response, varargin)
            %   stimulus    - vector of stimulus values
            %   response    - vector of responses
            
            if nargin < 2
                % simulating stimulus and response:
                warning('no stimulus & response provided-- SIMULATING!')
                warning('simulating a subject that has early weighting');
                nTrials             = 1e3;
                stimulus            = randn(nTrials, 7);
                subjectWeighting    = linspace(1.5, 0.5, 7);
                response            = rand(nTrials,1) < funpmf('logistic', mean(bsxfun(@times, stimulus, subjectWeighting),2), [0 10]);
            end
            
            p = inputParser();
            p.addOptional('useStimulusMean', false);
            p.addOptional('minimumTrialNumber', 10);
            p.addOptional('ridge', false);
            p.parse(varargin{:})
            
            assert(size(stimulus,1) == size(response(:),1), 'ERROR: number of rows in ''stimulus'' does not match number of responses');
            
            % get rid of nan trials:
            idxNan1             = isnan(response);
            idxNan2             = any(isnan(stimulus),2);
            self.response       = response(~(idxNan1 | idxNan2));
            self.stimulus       = stimulus(~(idxNan1 | idxNan2), :);
            self.nTrials        = size(self.response(:), 1);
            
            if self.nTrials < p.Results.minimumTrialNumber
                warning(['The number of trials in your dataset are smaller than the defined minimum of ' num2str(p.Results.minimumTrialNumber) ' trials. Aborting ppk...'])
                return;
            end
            nWeights            = size(self.stimulus,2);
            % design matrix along with indices for stimulus & bias:
            self.designMatrix   = [self.stimulus, ones(self.nTrials,1)];  % augmented by bias column
            % abbreviate and perform regression:
            X = self.designMatrix;
            Y = self.response;
            
            if p.Results.ridge
                rhovals=[0 .1 .2 .9 .99];
                rhoNull=.1;
                [wRidge, rho, w_se]  = autoRegress_logisticRidge(X, Y, 1:7, rhoNull, rhovals);
                self.costFunction = [];
                self.w = wRidge(1:nWeights);
                self.w_se = w_se(1:nWeights);
            else
                w0 = (X'*X)\(X'*Y);     % get least sqaured error weights as initial weights for minimization
                
                
                if any(isnan(w0))
                    warning(['you have NANs in your weight vector. big no no. Aborting ppk...'])
                    return;
                end
                
                opts                = optimset('GradObj', 'on', 'Hessian', 'on'); %#ok<*NASGU> % optimization options
                self.costFunction   = @(w) neglogli_bernoulliGLM(w, X, Y);   % negative log likelihood (requires nccLabCode)
                evalc('[w_mle] = fminunc(self.costFunction, w0, opts);');
                Hess                = compHess(self.costFunction, w_mle, 0.001);
                w_se          = sqrt(diag(inv(Hess)));
                % store weights+se for stimulus indices only, no bias:
                self.w        = w_mle(1:nWeights);
                self.w_se     = w_se(1:nWeights);
            end
            
            % compute l2 norm as well
            self.w_norm       = self.w ./ norm(self.w);
            self.w_se_norm    = self.w_se  ./ norm(self.w);
            
            % weight index - this is our standard weight index, designed
            % to work with 7 pulses (but omits pulse #4):
            % !!!!!
            % the index is evil!
            % (a) it omits pulse 4.
            % (b) it's not a real index!! because weight values can be
            % negaitve, the index is not always between -1 and 1, which
            % kinda defies the purpose of an index. so fuck it.
            if size(self.stimulus, 2)==7
                self.w_index      = (sum(self.w_norm(1:3)) - sum(self.w_norm(5:7))) / (sum(self.w_norm(1:3)) + sum(self.w_norm(5:7)));
            else
                self.w_index = nan;
            end
            
            % fit a linear model to the weights:
            self.w_lm = fitlm(1:numel(self.w_norm), self.w_norm);
            
            % compute center of mass (com) to weights:
            % !!!
            % this was the original method. it's wrong.
            %               self.w_com = median((self.w ./ sum(self.w)) .* (1:nWeights)');
            % the following method is a poorly coded (but works!):
            W = self.w_norm;
            X       = linspace(1, nWeights, nWeights);
            interpFactor = 1000;
            Xq = linspace(1, nWeights, nWeights * interpFactor);
            Wq = interp1(X, W, Xq);
            
            self.w_com = sum(Xq.*Wq)/sum(Wq);
            
            % weight_energy is the sum of squared residuals from the mean
            % of the weights. I perform this on the normalized weights
            self.w_energy = sum((self.w_norm - mean(self.w_norm)).^2);
        end
        
        
        
        %%
        function self = plot(self, varargin)
            %   self = plot(self, varargin)
            %
            % plots the weights over stimulus columns
            %
            % varargin -
            %   'plotFit' - set to true to plot the linear model fit to the
            %   weights (default is true).
            %   'normWeights' - set to true to plot the weights normalized
            %   by the l2 norm (default is true)
            
            if isempty(self.w)
                error('ERROR. weights have yet to be computed. cannot plot data');
            end
            
            figure(gcf); gca; hold on
            fh=plot(nan,nan);
            clr=get(fh, 'Color');   % get default color
            delete(fh);
            
            p = inputParser();
            p.addOptional('normWeights', true);
            p.addOptional('plotFit', true);
            p.addOptional('Color', clr);
            p.addOptional('xLabel', 'Pulse Number');
            p.addOptional('yLabel', 'Weight (a.u.)');
            p.addOptional('fontSize', 10)
            p.addOptional('Marker', '.')
            p.addOptional('MarkerSize', 3)
            p.addOptional('MarkerFaceColor', clr);
            p.addOptional('LineWidth', 1)
            p.addOptional('LineStyle', '-')
            p.addOptional('verbose', false)
            p.parse(varargin{:})
            
            clr          = p.Results.Color;
            marker       = p.Results.Marker;
            markerSize   = p.Results.MarkerSize;
            lineWidth    = p.Results.LineWidth;
            lineStyle    = p.Results.LineStyle;
            fontSize     = p.Results.fontSize;
            
            % x values for plotting.
            xValues = 1:size(self.stimulus,2);
            xlim([xValues(1)-0.5 xValues(end)+0.5])
            line(xlim, [0 0], 'LineStyle', '--', 'Color', [.5 .5 .5])
            
            % plot weights:
            if p.Results.normWeights
                hD = errorbarFancy(xValues, self.w_norm, self.w_se_norm, 'Color', clr, 'LineStyle', lineStyle, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
                ylim([-0.2 1])
            else
                hD = errorbarFancy(xValues, self.w, self.w_se, 'Color', clr, 'LineStyle', lineStyle, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
                ylim([-0.5 2])
            end
            
            % plot the linear fit:
            if p.Results.plotFit
                a       = table2array(self.w_lm.Coefficients(2,1));
                b       = table2array(self.w_lm.Coefficients(1,1));
                xFit    = [xValues(1)-0.5, xValues(end)+0.5];
                yFit    = a.*xFit + b;
                hF = plot(xFit, yFit, 'Color',  clr);
                
                self.hF = hF;
                
            end
            
            xticks([0 1 2 3 4 5 6 7 8])
            xlabel(p.Results.xLabel, 'fontSize', fontSize)
            ylabel(p.Results.yLabel, 'fontSize', fontSize)
            
            % add text to the plot:
            if p.Results.verbose
                text(min(xlim)+abs(.1*min(xlim)), min(ylim) + abs(.1*(min(ylim))), ['nTrials = ' num2str(self.nTrials)], 'Color', clr, 'fontsize', 8)
            end
            
            self.hD = hD;
            %self.hF = hF;
            
        end
        
        
        %%
        function self = plotDiagnostics(self, varargin)
            %   self = plotDiagnostics(self, varargin)
            %
            % plots a couple of useful things for diagnosing a problem
            
            if isempty(self.w)
                error('ERROR. weights have yet to be computed. cannot plot data');
            end
            
            figure; clf; hold on
            fh=plot(nan,nan);
            clr=get(fh, 'Color');   % get default color
            delete(fh);
            
            subplot(1,3,1)
            self.plot('verbose', true);
            hold on
            title('weights')
            
            subplot(1,3,2)
            nWeights = numel(self.w);
            imagesc(self.designMatrix(:, 1:nWeights))
            hold on
            xlabel('weight #')
            ylabel('trial #')
            title('design matrix')
            colorbar
            
            subplot(1,12,9)
            imagesc(mean(self.stimulus,2))
            hold on
            ylabel('trial #')
            title('stimulus mean')
            
            subplot(1,12,12)
            imagesc(self.response)
            hold on
            ylabel('trial #')
            title('responses')
            
        end
        
        %%
        function self = xValidate(self, varargin)
            %   self = xValidate(self, varargin)
            
            p = inputParser();
            p.addOptional('kFolds', 8);
            p.parse(varargin{:})
            
            % JAKE BUILD PRITTY PLEASE!
        end
        
    end
    
    methods(Access = private)
        
    end
    
end



%% HELPER FUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%

%%
function h = errorbarFancy(x,y,l,varargin)
%   h = errorbarFancy(x,y,l,u,varargin)
%
%   Error bar plot without the annoying horizontal lines.
%   Based on drawing a line, and then plotting.
%   takes in same arguments as the plot function, see LineSpec.

if nargin > 3 && isnumeric(varargin{1})
    u=varargin{1};
    if numel(varargin)>1
        varargin = varargin(2:end);
    else
        varargin={};
    end
else
    u=l;
end

% if ~exist('u', 'var') && exist('l', 'var')
%     u = l; % symmetrical errorbars
% end

% hack to find a lineWidth, if inputted:
if any(strcmp(varargin,'LineWidth'))
    idx = find(strcmp(varargin,'LineWidth'));
    lw  = varargin{idx+1};
else
    lw = 1; % default;
end

% hack to find a Color, if inputted:
if any(strcmp(varargin,'Color'))
    idx = find(strcmp(varargin,'Color'));
    clr  = varargin{idx+1};
else
    fh=plot(nan,nan);
    clr=get(fh, 'Color');
    delete(fh);
    varargin={varargin{:}, 'Color', clr};
    %     clr = [0 0 0]; % default;
end

gca; hold on;
line([x(:) x(:)]', [y(:)-l(:) y(:)+u(:)]', 'Color', clr, 'LineWidth', lw);
h = plot(x,y, 'Marker','o', 'LineStyle', 'None', varargin{:});


end

%%
function [H, g] = compHess(fun, x0, dx, varargin)
% [H, g] = compHess(fun, x0, dx, varargin)
% Numerically computes the Hessian of a function fun around point x0
% expects fun to have sytax:  y = fun(x, varargin);
%
% Input:
%   fun: @(x) function handle of a real valued function that takes column vector
%   x0: (n x 1) point at which Hessian and gradient are estimated
%   dx: (1) or (n x 1) step size for finite difference
%   extra arguments are passed to the fun
%
% Output:
%   H: Hessian estimate
%   g: gradient estiamte

n = numel(x0);
H = zeros(n,n);
g = zeros(n,1);
f0 = feval(fun, x0, varargin{:});

% input check
if ~all(isfloat(dx) & isfinite(dx))
    error('dx must be finite and float');
end
if any(dx <= 0)
    error('dx must be strictly positive');
end

if isscalar(dx)
    vdx = dx*ones(n,1);
elseif numel(dx) == n
    vdx = dx(:);
else
    error('vector dx must be the same size as x0');
end
A = diag(vdx/2);

for j = 1:n  % compute diagonal terms
    %central differences
    f1 = feval(fun, x0+2*A(:,j),varargin{:});
    f2 = feval(fun, x0-2*A(:,j),varargin{:});
    H(j,j) = f1+f2-2*f0;
    g(j) = (f1-f2)/2;
    % forward differences
    %     f1 = feval(fun, x0+2*A(:,j),varargin{:});
    %     f2 = feval(fun, x0+A(:,j),varargin{:});
    %     fx = feval(fun, x0,varargin{:});
    %     H(j,j) = f1-2*f2+fx;
    %     g(j)   = f2-fx;
end

for j = 1:n-1       % compute cross terms
    for i = j+1:n
        %central differences
        f11 = feval(fun, x0+A(:,j)+A(:,i),varargin{:});
        f22 = feval(fun, x0-A(:,j)-A(:,i),varargin{:});
        f12 = feval(fun, x0+A(:,j)-A(:,i),varargin{:});
        f21 = feval(fun, x0-A(:,j)+A(:,i),varargin{:});
        H(j,i) = f11+f22-f12-f21;
        % forward differences
        %         fx = feval(fun, x0,varargin{:});
        %         f1 = feval(fun, x0+A(:,j)+A(:,i),varargin{:});
        %         f12 = feval(fun, x0+A(:,j),varargin{:});
        %         f21 = feval(fun, x0+A(:,i),varargin{:});
        %         H(j,i) = f1-f12-f21-fx;
        
        H(i,j) = H(j,i);
    end
end

H = H./(vdx * vdx');
g = g./vdx;

end

%%
function [L,dL,ddL] = neglogli_bernoulliGLM(wts,X,Y)
% [L,dL,ddL] = neglogli_bernoulliGLM(wts,X,Y)
%
% Compute negative log-likelihood of data under logistic regression model,
% plus gradient and Hessian
%
% Inputs:
% wts [m x 1] - regression weights
%   X [N x m] - regressors
%   Y [N x 1] - output (binary vector of 1s and 0s).

xproj = X*wts;

if nargout <= 1
    L = -Y'*xproj + sum(softrect(xproj)); % neg log-likelihood
    
elseif nargout == 2
    [f,df] = softrect(xproj); % evaluate log-normalizer
    L = -Y'*xproj + sum(f); % neg log-likelihood
    dL = X'*(df-Y);         % gradient
    
elseif nargout == 3
    [f,df,ddf] = softrect(xproj); % evaluate log-normalizer
    L = -Y'*xproj + sum(f); % neg log-likelihood
    dL = X'*(df-Y);         % gradient
    ddL = X'*bsxfun(@times,X,ddf); % Hessian
end


end

%%
% -------------------------------------------------------------
% ----- SoftRect Function (log-normalizer) --------------------
% -------------------------------------------------------------

function [f,df,ddf] = softrect(x)
%  [f,df,ddf] = softrect(x);
%
%  Computes: f(x) = log(1+exp(x))
%  and first and second derivatives

f = log(1+exp(x));

if nargout > 1
    df = exp(x)./(1+exp(x));
end
if nargout > 2
    ddf = exp(x)./(1+exp(x)).^2;
end

% Check for small values
if any(x(:)<-20)
    iix = (x(:)<-20);
    f(iix) = exp(x(iix));
    df(iix) = f(iix);
    ddf(iix) = f(iix);
end

% Check for large values
if any(x(:)>500)
    iix = (x(:)>500);
    f(iix) = x(iix);
    df(iix) = 1;
    ddf(iix) = 0;
end

end
%%

function [wRidge,rho,SDebars,postHess,logevid] = autoRegress_logisticRidge(xx,yy,iirdge,rhoNull,rhovals,w0)
% [wRidge,rho,SDebars,postHess,logevid] = autoRegress_logisticRidge(xx,yy,iirdge,rhoNull,rhovals,w0)
%
% Computes empirical Bayes logistic ridge regression filter estimate under a
% Bernoulli-GLM with ridge prior
%
% Inputs:
%        xx [n x p] - stimulus (regressors)
%        yy [n x 1] - spike counts (response)
%     iirdg [q x 1] - weight indices for which to learn ridge parameter (columns of xx)
%   rhoNull [1 x 1] - fixed prior precision for other columns of xx
%   rhovals [v x 1] - list of prior precisions to use for grid search
%        w0 [p x 1] - initial estimate of regression weights (optional)
%   opts (optional) = options stucture: fields 'tolX', 'tolFun', 'maxIter','verbose'
%
% Outputs:
%   wRidge [p x 1] - empirical bayes estimate of weight vector
%      rho [1 x 1] - estimate for prior precision (inverse variance)
%  SDebars [p x 1] - 1 SD error bars from posterior
% postHess [p x p] - Hessian (2nd derivs) of negative log-posterior at
%                    maximum ( inverse of posterior covariance)
%  logevid [1 x 1] - log-evidence for hyperparameters
%
% $Id: autoRegress_logisticRidge.m 2074 2012-09-18 16:55:46Z pillow $

nw = size(xx,2); % length of stimulus filter

% initialize filter estimate with MAP regression estimate, if necessary
if nargin == 5
    rho0 = 5;        % initial value of ridge parameter
    Lprior = speye(nw)*rho0;
    w0 = (xx'*xx+Lprior)\(xx'*yy);
end

% --- set prior and log-likelihood function pointers ---
mstruct.neglogli = @neglogli_bernoulliGLM;
mstruct.logprior = @logprior_ridge;
mstruct.liargs = {xx,yy};
mstruct.priargs = {iirdge,rhoNull};

% --- Do grid search over ridge parameter -----
[hprsMax,wmapMax] = gridsearch_GLMevidence(w0,mstruct,rhovals);
fprintf('best grid point: rho (precision)=%.1f\n', hprsMax);

% --- Do gradient ascent on evidence ----
[wRidge,rho,logevid,postHess] = findEBestimate_GLM(wmapMax,hprsMax,mstruct);

if nargout >2
    SDebars = sqrt(diag(inv(postHess)));
end

end

function [wEB,hprsML,logEv,postHess] = findEBestimate_GLM(w0,hprs0,mstruct)
% [wEB,hprsML,logEv,postHess] = findEBestimate_GLM(w0,hprs0,mstruct)
%
% Find empirical Bayes estimate for GLM regression weights (Poisson or
% logistic). Maximizes marginal likelihood (computed using Laplace
% approximation) for hyperparameters and returns MAP estimate given those
% hyperparameters
%
% INPUTS:
%      w0 [m x 1] - initial guess at parameters
%   hprs0 [p x 1] - initial guess at hyper-parameters
% mstruct [1 x 1] - model structure with fields:
%         .neglogli - func handle for negative log-likelihood
%         .logprior - func handle for log-prior
%         .liargs - cell array with args to neg log-likelihood
%         .priargs - cell array with args to log-prior function
%
% OUTPUTS:
%        w [m x 1] - empirical Bayes parameter estimate for regression weights
%     hprs [p x 1] - maximum marginal likelihood hyper-parameter estimate
%    logEv [1 x 1] - log-evidence at maximum
% postHess [1 x 1] - Hessian (2nd derivs) of negative log-posterior at
%                    maximum ( inverse of posterior covariance)
%
% Note: this implementation is not especially robust or efficient, and
% should only be used when initialized somewhere in the vicinity of the ML
% hyperparameter values (e.g., identified by a coarse grid search).
%
% Parametrization: log transforms first hyperparameter (assumed to be precision),
% and logit transforms remaining hyperprs (assumed to be in [0,1]).
%
% $Id: findEBestimate_GLM.m 2052 2012-09-11 22:08:48Z pillow $

% find MAP estimate for w given initial hyperparams
opts1 = struct('tolX',1e-10,'tolFun',1e-10,'maxIter',1e4,'verbose',0);
lfpost = @(w)(neglogpost_GLM(w,hprs0,mstruct)); % loss func handle
wmap0 = fminNewton(lfpost,w0,opts1); % do optimization for w_map

% transform hyperparameters to Reals
htprs0 = transformhprs(hprs0,1);

% Find maximum of evidence as a function of hyperparams
opts2 = struct('tolX',1e-10,'tolFun',1e-10,'maxIter',25,'verbose',0); % opts governing MAP search
hloss = @(ht)(updateMAPandEvalEvidence(ht,wmap0,opts2,mstruct)); % neg evidence func handle
opts_fminunc = optimset('display', 'iter', 'largescale', 'off','maxfunevals',1e3,...
    'FinDiffType','central'); % opts governing evidence search
[htprsML,neglogEv1] = fminunc(hloss,htprs0,opts_fminunc); % maximize evidence

% get MAP estimate at these hyperparams
[neglogEv2,wEB,postHess] = updateMAPandEvalEvidence(htprsML,wmap0,opts1,mstruct);
if abs(neglogEv1-neglogEv2)>.1
    warning('evidence estimate may be inaccurate (findEBestimate_bernoulliGLM)');
end

% un-transform hyperparameters to appropriate range
hprsML = transformhprs(htprsML,-1);
logEv = -neglogEv2;

end

% ======================================================================
% Sub-function to update MAP parameters and evaluate evidence
% ======================================================================

function [neglogEv,wmap,postHess] = updateMAPandEvalEvidence(htprs,wmap0,opts,mstruct)
% [neglogEv,wmap] = updateMAPandEvalEvidence(htprs,X,Y,fprior,wmap0,opts)

% un-transform hyperparameters
hprs = transformhprs(htprs,-1);

% find MAP estimate for w
lfpost = @(w)(neglogpost_GLM(w,hprs,mstruct));
wmap = fminNewton(lfpost,wmap0,opts);

% evaluate negative log evidence using Laplace approxmiation
neglogEv = -logevid_GLM(wmap,hprs,mstruct);

% Compute Hessian of negative log-posterior
if nargout > 2
    [~,~,postHess] = lfpost(wmap);
end

end