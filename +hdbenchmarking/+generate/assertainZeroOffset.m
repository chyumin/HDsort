function [fp_out, P] = assertainZeroOffset(fp_in, varargin)

P.tukey_ratio = 0.1;
P.subtractMean = true;
P = util.parseInputs(P, varargin, 'error');

[nTf, nCh, nU] = size(fp_in);

P.meanPerChannelPerUnit = mean(fp_in, 1);

if P.subtractMean
    fp_ = fp_in - repmat(P.meanPerChannelPerUnit, nTf, 1, 1);
else
    fp_ = fp_in;
    
    if any(abs(P.meanPerChannelPerUnit) > 0.1)
        warning('The input waveform contains a DC component')
    end
end

twindow_ = tukeywin(nTf, P.tukey_ratio);
P.twindow = repmat(twindow_, 1, nCh, nU);

fp_out = fp_ .* P.twindow;
