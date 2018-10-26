

ST = {[100, 120, 140]',
    [112, 229, 400, 410]',
    [95, 290]',
    [10, 120, 121, 122, 411]'};

window = 0.001;
sampling = 20000;

[binnedST, edges, times] = hdsort.spiketrain.binSpikeTrains(ST, window, sampling);

%%
figure; hold on;
for ii = 1:4
    %subplot(2,2,ii)
    bar(times, binnedST(ii, :))
end





