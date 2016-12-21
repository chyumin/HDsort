
nC = 10;
X = randn(1000,nC);
X(2:end-1,2) = X(2:end-1,2) + X(1:end-2,1) + X(3:end, 3);

X(10:end, 7) = X(1:end-9, 4);
xc = xcorr(X, 10, 'none');
Cte = mysort.hdsort.noise.xcorr2Cte(xc);

Ccol = mysort.hdsort.noise.Cte2Ccol(Cte, nC);

xc2 = mysort.hdsort.noise.ccol2xcorr(Ccol);

Cte2 = mysort.hdsort.noise.ccol2Cte(Ccol);

figure;
subhdsort.plot.2,3,1);
imagesc(xc);
title('Matlab xcorr');

subhdsort.plot.2,3,4);
imagesc(Cte);
title('Time embedding');

subhdsort.plot.2,3,[2 5]);
imagesc(Ccol);
title('Channel embedding ccol');

subhdsort.plot.2,3,3);
imagesc(xc2);
title('Matlab xcorr 2');

subhdsort.plot.2,3,6);
imagesc(Cte2);
title('Time embedding 2');

norm(xc-xc2)
norm(Cte-Cte2)