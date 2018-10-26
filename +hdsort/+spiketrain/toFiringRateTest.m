%%
sampling = 0.01; 
window  = 0.1; 
timePeriod = [0, 1];
ST = [1000:10:1100, 2000:200:5000, 15000:1000:19000];
[fr, x] = hdsort.spiketrain.toFiringRate(ST/20000, window, sampling, timePeriod);
sum(fr)

figure
ah = subplot(2,1,1);
RP = hdsort.plot.Rasterplot({ST}, 'ah', ah)
ah(2) = subplot(2,1,2);
plot(x, fr)
linkaxes(ah, 'x');

