    nS = 100;
    nC = 2;

    X = repmat([5:-.5:0 1:10 11 10:-1:1 0:.5:5], nS, nC)-5;
    maxIdx = 22;
    X = X + 1*randn(size(X));
    shifts = [-4:2:4];
    stau = repmat(shifts', nS/length(shifts),1);
    figure; plot(X(1,:));
    X = hdsort.util.shiftMCRows(X,stau,nC,1);
    hold on; plot(X(1,:), 'r');
    debug = 1;

    
    [tau Y] = hdsort.waveforms.vAlignOnAverageMaxSample(X, nC, 'maxIdx', 22, 'nMaxChannelsForWeighting', 5);
    fprintf('Shift Error: %.4f\n', sum(abs(-tau-stau)));
    figure; plot(stau, -tau, '.');
    xlabel('real tau');
    ylabel('estimated tau');

    title('Aligned Spikes');    

