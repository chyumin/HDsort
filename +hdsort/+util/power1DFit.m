function [out] = power1DFit(X, Y, varargin)
% Input:
%   exp1DFit(x, y, varargin):
%  - X: datapoints on x-axis
%  - Y: datapoints on y-axis
%
%   exp1DFit(y, varargin):
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
P.y0Constraint = [];
P.bMaximum = [];
P.bMinimum = [];
P = hdsort.util.parseInputs(P, varargin, 'error');

% Create handle for parameter minimising function:
EH = @(params) E(params, X, Y);

% Estimate starting parameters:
%[y_min idxy_min] = min(Y);
%[y_max idxy_max] = max(Y);
[x_min idxx_min] = min(X);
[x_max idxx_max] = max(X);
if x_min < 0
    b_in = 1;
    a_in = sign(Y(idxx_max) - Y(idxx_min));
else
    b_in = -1;
    a_in = -sign(Y(idxx_max) - Y(idxx_min));
end
mu_in = 0;
y0_in = 0;

params0 = [a_in y0_in b_in mu_in];

% Run fit:
if ~isempty(P.MaxIter)
    [paramsopt, fopt, exitflag] = fminsearch(EH, params0, optimset('MaxIter', P.MaxIter));
else
    [paramsopt, fopt, exitflag] = fminsearch(EH, params0); 
end


    function [a, y0, b, mu] = interpretParameters(params)
        a = params(1);
        y0 = params(2);
        b = params(3);
        mu = params(4);
        
        %% Define Constraints:
        
        if ~isempty(P.muConstraint)
           mu = P.muConstraint;
        end
        if ~isempty(P.y0Constraint)
            y0 = P.y0Constraint;
        end
        
        if ~isempty(P.bMinimum)
            b = max(b, P.bMinimum);
        end
        if ~isempty(P.bMaximum)
            b = min(b, P.bMaximum);
        end
    end

    function e = E(params, X, Y)
        [a, y0, b, mu] = interpretParameters(params);
        yFit = yFitModel(X, a, y0, b, mu);
        e = norm(yFit - Y);
    end

    function yFit = yFitModel(X, a, y0, b, mu)
        yFit = a*(X-mu).^b + y0;
    end

[a_opt, y0_opt, b_opt, mu_opt] = interpretParameters(paramsopt);
out.yFit = yFitModel(X, a_opt, y0_opt, b_opt, mu_opt);

out.Yin = Y;
out.Xin = X;
out.errors = Y - out.yFit;

out.a = a_opt;
out.y0 = y0_opt;
out.b = b_opt;
out.mu = mu_opt;

out.fopt = fopt;
out.success = exitflag;
out.yFitModel = @(X) yFitModel(X, out.a, out.y0, out.b, out.mu);
out.inputParameters = P;

if P.debug
    figure; hold on;
    scatter(X, Y);
    
    x = sort(X);
    hdsort.plot.x, out.yFitModel(x), 'r');
    axis([0 1.1*max(X) 0 1.1*max(Y)])
end

end