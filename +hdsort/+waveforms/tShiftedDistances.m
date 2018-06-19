function [out, P] = tShiftedDistances(T1, T2, varargin)
% T1 and T2 are multi-electrode waveforms of several spikes or units.
% The output contains the distances between the pairs of units at each
% shift, as well as the minimal and maximal distances and their shift value
% (no shift: 0).

P.debug = false;
P.distance = 'euclidean';
P.maxShift = 10;
P = hdsort.util.parseInputs(P, varargin, 'error');

[nTf1, nCh1, nU1] = size(T1);
[nTf2, nCh2, nU2] = size(T2);

assert(nTf1 == nTf2 & nCh1 == nCh2, 'Number of channels and lengths of waveforms must be the same!')
%assert(ndims(T1) == 2 & ndims(T2) == 2, 'Input waveforms must be 2D!')

out.distance = zeros(P.maxShift, nU1, nU2);
%out.shifts = 0:P.maxShift;
out.shifts = -P.maxShift:P.maxShift;

for ii = 1:numel(out.shifts)
    tau = out.shifts(ii);
    T2_shifted = hdsort.waveforms.tShift(T2, tau, true);
    
    T1_ = reshape(T1, nCh1*nTf1, nU1)';
    T2_ = reshape(T2_shifted, nCh2*nTf2, nU2)';
    
    out.distance(ii, :, :) = pdist2( T1_, T2_, P.distance );
end

for ii = 1:nU1
    for jj = 1:nU2
        
        [out.minDistance(ii, jj), minShiftIdx] = min(out.distance(:, ii, jj));
        out.minShift(ii, jj) = out.shifts(minShiftIdx);
        
        [out.maxDistance(ii, jj), maxShiftIdx] = max(out.distance(:, ii, jj));
        out.maxShift(ii, jj) = out.shifts(maxShiftIdx);
    end
end

end