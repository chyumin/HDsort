function FTO = getFootprintOverlap(FP1, FP2)

if nargin == 1
    T = FP1;
else
    T = [];
    T(:,:,1) = FP1;
    T(:,:,2) = FP2;    
end
N = size(T, 3);
assert( N > 1, 'You can only compare more than one footprint!')
assert( N < 100, 'Otherwise it takes too long in order to calculate the XCorr!')

[xc A S] = hdsort.waveforms.tXCorr(T);

FTO = zeros(N,N);
for ii = 1:N
    for jj = 1:N

if S(jj,ii) > 0
    n = S(jj, ii);
    FP1_shifted = T(:,:,ii);
    FP2_shifted = [T((n+1):end, :, jj); zeros(n,size(T, 2))]; %[FP2((n+1):end, :); zeros(n,size(FP2, 2))];
else
    n = S(ii,jj);
    FP1_shifted = [T((n+1):end, :, ii); zeros(n,size(T, 2))];%[FP1((n+1):end, :); zeros(n,size(FP1, 2))];
    FP2_shifted = T(:,:,jj);
end

        FTO(ii, jj) = 100 * FP2_shifted(:)'*FP1_shifted(:) / (norm(FP2_shifted(:)) * norm(FP1_shifted(:)) );
    end
end

if nargin ~= 1
    assert(round(FTO(1,2)*10) == round(FTO(2,1)*10), 'Error here!')
    FTO = FTO(1,2);
end

if false
    %%
    figure; hold on;
    hdsort.plot.reshape(FP1, 1, []), 'b');
    hdsort.plot.reshape(FP2, 1, []), 'r');
    
    hdsort.plot.reshape( FP1_shifted, 1, [])+10, 'g');
    hdsort.plot.reshape( FP2_shifted, 1, [])+10, 'm');
end


