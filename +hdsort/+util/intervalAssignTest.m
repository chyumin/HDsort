
intervals = [1 4; 4 7; 7 12; 12 40]
hdsort.spiketrain = [0.2; 2.3; 11.9; 4.4; 10.3; 12; 10.4; 39.99; 40.0; 40.1]

[intervalIndices] = hdsort.util.intervalAssign(hdsort.spiketrain, intervals)
