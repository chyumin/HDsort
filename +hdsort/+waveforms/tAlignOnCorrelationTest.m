% Templates
T = [];
T(:,1,1) = [0 0  0 -1   -2 -1    4 5 4 3 2 1 0 0];
T(:,1,2) = T(:,1,1)/2;
T(:,1,3) = flipud(T(:,1,1));
T(:,1,4) = [0 0 -1 -1.5 -2 -1.5 .5 1 4 5 1 0 0 0];

xvsf = mysort.hdsort.util.calculateXIvsF(hdsort.waveforms.t2v(T),hdsort.waveforms.t2v(T),1,1);


[ali tau xvsfali] = hdsort.waveforms.tAlignOnCorrelation(T);
[ali2 tau2 xvsfali2] = hdsort.waveforms.tAlignOnCorrelation(ali);
figure;
subhdsort.plot.2,1,1)
hdsort.plot.squeeze(T));
subhdsort.plot.2,1,2)
hdsort.plot.squeeze(ali));


mysort.hdsort.plot.XIvsF(T, T)
mysort.hdsort.plot.XIvsF(ali, ali)
