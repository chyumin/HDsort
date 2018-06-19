function T = tShift(T, tau, trunc)
% shifts the single channel waveforms stored in T
% truncates in a way, that A has the same dimensions as T
% If a waveform has a shift of 0, and trunc is 1 it will be unchanged.

if nargin < 3; trunc = 0; end
tau = tau(:);
assert(~any(tau~=round(tau)), 'tau must contain only integers!');

[nTf, nC, nU] = size(T);

if size(tau, 1) == 1
    tau = repmat(tau', nC*nU, 1);
else
    assert(size(tau, 1) == nU, 'You must provide one tau for each unit!')
    tau = repmat(tau', nC, 1);
end
tau = tau(:);

T = hdsort.waveforms.t2m(T);
T = hdsort.util.shiftRows(T, tau, trunc);
T = hdsort.waveforms.m2t(T, nC);

end
