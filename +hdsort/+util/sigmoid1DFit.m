function [out] = sigmoid1DFit(X, Y, varargin)
% Input:
%   sigmoid1DFit(x, y, varargin):
%  - X: datapoints on x-axis
%  - Y: datapoints on y-axis
%
%   sigmoid1DFit(y, varargin):
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
P.yRangeConstraint = [];
P.y0Constraint = [];
P = hdsort.util.parseInputs(P, varargin, 'error');

% Create handle for parameter minimising function:
EH = @(params) E(params, X, Y);

% Estimate starting parameters:
xRange_in = max(X) - min(X);
yRange_in = max(Y) - min(Y);
y0_in = min(Y);
sigma_in = yRange_in / xRange_in;
mu_in = mean(X);
params0 = [yRange_in y0_in sigma_in mu_in];

% Run fit:
if ~isempty(P.MaxIter)
    [paramsopt, fopt, exitflag] = fminsearch(EH, params0, optimset('MaxIter', P.MaxIter));
else
    [paramsopt, fopt, exitflag] = fminsearch(EH, params0); 
end


    function [yRange, y0, sigma, mu] = interpretParameters(params)
        yRange = params(1);
        y0 = params(2);
        sigma = params(3);
        mu = params(4);
        
        %% Define Constraints:
        
        if ~isempty(P.yRangeConstraint)
            yRange = P.yRangeConstraint;
        end
        if ~isempty(P.y0Constraint)
            y0 = P.y0Constraint;
        end
        
    end

    function e = E(params, X, Y)
        [yRange, y0, sigma, mu] = interpretParameters(params);
        yFit = yFitModel(X, yRange, y0, sigma, mu);
        e = norm(yFit - Y);
    end

    function yFit = yFitModel(X, yRange, y0, sigma, mu)
        yFit = yRange ./ (1.0 + exp( -sigma * (X - mu) )) + y0;     
    end

[yRange_opt, y0_opt, sigma_opt, mu_opt] = interpretParameters(paramsopt);
out.yFit = yFitModel(X, yRange_opt, y0_opt, sigma_opt, mu_opt);

out.Yin = Y;
out.Xin = X;
out.errors = Y - out.yFit;

out.yRange = yRange_opt;
out.y0 = y0_opt;
out.sigma = sigma_opt;
out.mu = mu_opt;

out.fopt = fopt;
out.success = exitflag;
out.yFitModel = @(X) yFitModel(X, out.yRange, out.y0, out.sigma, out.mu);
out.inputParameters = P;

if P.debug
    figure; hold on;
    scatter(X, Y);
    
    x = min(X):max(X);
    hdsort.plot.x, out.yFitModel(x), 'r');
end

end