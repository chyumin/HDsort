function CCol = xcorr2CCol(xc)
%     maxlag = (size(xc,1)-1)/2;
%     Tf = maxlag+1;
    nC = sqrt(size(xc,2));
    Cte = hdsort.noise.xcorr2Cte(xc);
    CCol = hdsort.noise.Cte2Ccol(Cte, nC);