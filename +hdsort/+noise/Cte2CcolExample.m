
nC = 11;
X = randn(1000,nC);
X(2:end-1,2) = X(2:end-1,2) + X(1:end-2,1) + X(3:end, 3);

X(10:end, 7) = X(1:end-9, 4);
xc = xcorr(X, 10, 'none');
Cte = mysort.hdsort.noise.xcorr2Cte(xc);

Ccol = mysort.hdsort.noise.Cte2Ccol(Cte, nC);

figure;
subhdsort.plot.2,2,1);
imagesc(xc);
title('Matlab xcorr');

subhdsort.plot.2,2,3);
imagesc(Cte);
title('Time embedding');

subhdsort.plot.2,2,[2 4]);
imagesc(Ccol);
title('Channel embedding ccol');