function [out] = gaussian1DFit(X, Y, varargin)
% Input:
%   gaussian1DFit(x, y, varargin):
%  - X: datapoints on x-axis
%  - Y: datapoints on y-axis
%
%   gaussian1DFit(y, varargin):
%  - XY: two column matrix with datapoints in (x,y)-format
%
%   varargin: 
%  - 'MaxIter': maximal number of iterations
%
% Output:
%  - out: structure containing all fitting parameters

if nargin == 1
    assert(size(X,2) == 2, 'If only one input is privided, it must be a 2-column vector!')
    Y = X(:,2);
    X = X(:,1);
    varargin = {};
elseif ischar(Y)
    assert(size(X,2) == 2, 'If only one input is privided, it must be a 2-column vector!')
    Y = X(:,2);
    X = X(:,1);
    varargin = {Y, varargin{:}};
end
assert(size(X,1) == size(Y,1), 'X and Y must have the same length!')

P.debug = false;
P.MaxIter = [];
P.muConstraint = [];
P.sigmaConstraint = [];
P.ampConstraint = [];
P = hdsort.util.parseInputs(P, varargin, 'error');

% Create handle for parameter minimising function:
EH = @(params) E(params, X, Y);

% Estimate starting parameters:
amp_in = max(X);
sigma_in = std(X);
mu_in = mean(X);

params0 = [amp_in, sigma_in, mu_in];

% Run fit:
if ~isempty(P.MaxIter)
    [paramsopt, fopt, exitflag] = fminsearch(EH, params0, optimset('MaxIter', P.MaxIter));
else
    [paramsopt, fopt, exitflag] = fminsearch(EH, params0); 
end

    function [amp, sigma, mu] = interpretParameters(params)
        amp = params(1);
        sigma = params(2);
        mu = params(3);
        
        %% Define Constraints:
        if ~isempty(P.ampConstraint)
            y0 = P.ampConstraint;
        end
        if ~isempty(P.muConstraint)
           mu = P.muConstraint;
        end
        if ~isempty(P.sigmaConstraint)
           mu = P.sigmaConstraint;
        end
    end
    
    function e = E(params, X, Y)
        [amp, sigma, mu] = interpretParameters(params);
        yFit = yFitModel(X, amp, sigma, mu);
        e = norm(yFit - Y);
    end

    function yFit = yFitModel(X, amp, sigma, mu)
        yFit = amp * gaussmf(X, [sigma, mu]);
    end

[amp_opt, sigma_opt, mu_opt] = interpretParameters(paramsopt);
out.yFit = yFitModel(X, amp_opt, sigma_opt, mu_opt);

out.Yin = Y;
out.Xin = X;
out.errors = Y - out.yFit;

out.amp = amp_opt;
out.sigma = sigma_opt;
out.mu = mu_opt;

out.fopt = fopt;
out.success = exitflag;
out.yFitModel = @(X) yFitModel(X, out.amp, out.sigma, out.mu);
out.inputParameters = P;



if P.debug
    %%
    figure; hold on;
    scatter(X, Y);
    
    X_ = linspace( out.mu-5*out.sigma, out.mu+5*out.sigma, 100);
    plot(X_, out.yFitModel(X_), 'r');
    axis([1.1*min(X_) 1.1*max(X_), 0, 1.1*out.amp]);
end

end