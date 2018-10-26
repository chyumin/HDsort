% This tutorial explains general functions that can be used to access data
% from files and to perform very simple data analysis.
% 
% For a tutorial on how to use the spike-sorting functions, please read
% 'sorting_example.m'

%% Add 'External' functions to path:
setup()

%% Download test file:
mainFolder = fullfile('.', 'Tutorial');
fileName = 'chirpsweep01.raw.h5';
rawFile = fullfile(mainFolder, fileName);
download_testfile(rawFile)
assert(exist(rawFile, 'file')>0, 'This tutorial requires a raw data file!')

%% Get a handle onto the raw file, which will then be used as input to the sorter:
sortingName = 'sorting_example';
RAW = hdsort.file.BELMEAFile(rawFile);
RAW.restrictToConnectedChannels(); % This line is very important to not sort empty electrodes!

%% Extract data from the file:
sampleIdx = 1:20000;
channelIdx = 1:10;

Xunfiltered = RAW.getData(sampleIdx, channelIdx);
Xfiltered = RAW.getFilteredData(sampleIdx, channelIdx);

figure; hold on
plot(Xfiltered);
plot(Xunfiltered);

%% Plot the electrode configuration with which the data was recorded:
RAW.MultiElectrode.plotConfig();

%% Load data into memory, filter it and create a wrapper around it to use advanced functions:
sampleIdx = 1:(20000*60);
channelIdx = 1:10;

RAW.restrictToChannels(channelIdx); 
X = RAW.getFilteredData(sampleIdx, :);
M = hdsort.file.DataMatrix(X, 'MultiElectrode', RAW.MultiElectrode);

%% Detect spikes:
spikeTrains = M.detectSpikes('thr', 5);

%% Cut all spikes that were detected in one channel:
channel = 7;
spiketimes = spikeTrains{channel};
wf = M.getWaveform(spiketimes, 20, 50);
twf = hdsort.waveforms.v2t(wf, size(M, 2)); % reshape the data

%% Plot an average waveform "footprint":
median_wf = median(twf, 3);
hdsort.plot.Waveforms2D(median_wf, M.MultiElectrode.electrodePositions);

