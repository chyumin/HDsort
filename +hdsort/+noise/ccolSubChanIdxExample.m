
nC = 10;
X = randn(1000,nC);
X(2:end-1,2) = X(2:end-1,2) + X(1:end-2,1) + X(3:end, 3);

X(10:end, 7) = X(1:end-9, 4);
xc = xcorr(X, 10, 'none');
Cte = mysort.hdsort.noise.xcorr2Cte(xc);

Ccol = mysort.hdsort.noise.Cte2Ccol(Cte, nC);

Ccol_red = mysort.hdsort.noise.ccolSubChanIdx(Ccol, [2:3 4 7]);

xc_red = mysort.hdsort.noise.ccol2xcorr(Ccol_red);

Cte_red = mysort.hdsort.noise.ccol2Cte(Ccol_red);

figure;
subhdsort.plot.2,4,1);
imagesc(xc);
title('Matlab xcorr');

subhdsort.plot.2,4,5);
imagesc(Cte);
title('Time embedding');

subhdsort.plot.2,4,[2 6]);
imagesc(Ccol);
title('Channel embedding ccol');

subhdsort.plot.2,4,[3 7]);
imagesc(Ccol_red);
title('Channel embedding ccol');

subhdsort.plot.2,4,4);
imagesc(xc_red);
title('Matlab xcorr 2');

subhdsort.plot.2,4,8);
imagesc(Cte_red);
title('Time embedding 2');
