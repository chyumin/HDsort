function [X tau] = vAlignOnTemplateInterpolated(X,nC,t,maxShift)
    if nargin < 4
        maxShift = 5;
    end
    opt = [];
    opt.MaxFunEvals = 30;
    opt.TolX = 0.01;
    opt.Display = 'off';    
    
    nS = size(X,1);

    mt = hdsort.waveforms.v2m(t, nC);
    range = 1:size(mt,2);
    
    tau = zeros(nS, 1);
    for i=1:nS
        xsinc = hdsort.waveforms.mSincfun(hdsort.waveforms.v2m(X(i,:), nC));            
        fun = @(t)   -sum(sum(mt.*xsinc(range-t)));
        [tau(i), k] = fminbnd(fun, -maxShift, maxShift, opt);
    end

    X   = hdsort.util.shiftRowsInterpolated(X, tau, nC);
end