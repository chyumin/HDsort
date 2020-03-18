function P = defaultParameters()

P = struct();

P.spikeDetection.method = '-';
P.spikeDetection.thr = 4.2;
P.artefactDetection.use = 0;

P.legs.maxElPerGroup = 9;
P.legs.minElPerGroup = 1;
P.legs.addIfNearerThan = 20; % always add direct neighbors
P.legs.maxDistanceWithinGroup = 52;

P.botm.run = 0;
P.botm.Tf = 75;
P.botm.cutLeft = 20;

P.spikeCutting.maxSpikes = 200000000000; % Set this to basically inf
P.spikeCutting.blockwise = false;

P.noiseEstimation.minDistFromSpikes = 80;

P.spikeAlignment.initAlignment = '-';
P.spikeAlignment.maxSpikes = 50000;     % so many spikes will be clustered

P.featureExtraction.nDims = 6;

P.clustering.maxSpikes = 50000;  % dont align spikes you dont cluster...
P.clustering.meanShiftBandWidthFactor = 1.8;
%P.clustering.meanShiftBandWidth = sqrt(1.8*6); % todo: check this!

P.mergeTemplates.merge = 1;
P.mergeTemplates.upsampleFactor = 3;
P.mergeTemplates.atCorrelation = .93; % DONT SET THIS TOO LOW! USE OTHER ELECTRODES ON FULL FOOTPRINT TO MERGE
P.mergeTemplates.ifMaxRelDistSmallerPercent = 30;


P.templateEstimation.cutLeft = 10;
P.templateEstimation.Tf = 55;
P.templateEstimation.maxSpikes = 100;

end

% %% hdsort.sortScript.m Parameters:
% P.artefactDetection.use = 0;
% P.artefactDetection.width = 350;
% P.artefactDetection.threshold = size(DS,2)*41; % should be sampling frequency dependent!
% 
% % Spike Detection
% P.spikeDetection.method = '-';
% P.spikeDetection.thr = 3.5;
% P.spikeDetection.minDist = 25;
% P.spikeDetection.maxDataLength = [];
% P.spikeDetection.mergeSpikesMaxDist = 18; % do not go below 13 for HDMEA!
% P.spikeDetection.removeEventsWithAbsAmplitudeLargerThan = 1500000;
% P.spikeDetection.removeOscillationGroupsWithNMembersOfMinAmplitude = [0 0 0]; % should remove atrifacts that are oscillations (successive detections of peaks with alternating signs). First parameter is N, second the minHeight, third the maxDist within a group
% 
% P.spikeCutting.maxSpikes = 200000000000;
% 
% P.spikeCutting.Tf = 71;
% P.spikeCutting.cutLeft = 20;
% P.spikeCutting.chunkSize = 10000;
% P.spikeCutting.blockwise = true;
% 
% % Spike Alignment
% P.spikeAlignment.method = 'onUpsampledMean';
% P.spikeAlignment.maxSpikes = 50000;
% P.spikeAlignment.Tf = 35;
% P.spikeAlignment.cutLeft = 12;
% P.spikeAlignment.initAlignment = []; %'-' % Do not use init alignment
% P.spikeAlignment.maxIdx = P.spikeAlignment.cutLeft + 1;
% 
% P.spikeAlignment.maxIterations = 30;
% 
% % Noise estimation
% P.noiseEstimation.maxSamples = 1000000;
% P.noiseEstimation.maxSamplesPerEpoch = 500000;
% P.noiseEstimation.minDistFromSpikes = 60;
% 
% % Feature extraction
% P.featureExtraction.Tf = 20;
% P.featureExtraction.cutLeft = 8;
% P.featureExtraction.nDims = 6;
% 
% % Clustering
% P.clustering.maxSpikes = 30000;
% P.clustering.meanShiftBandWidth = []; % default is sqrt( P.clustering.meanShiftBandWidth * P.featureExtraction.nDims )
% P.clustering.meanShiftBandWidthFactor = 1.1;
% P.clustering.minSpikesPerCluster = 10;
% 
% % Template Matching on Cut Spikes
% P.templateMatchingCut.prior = .0001;
% P.templateMatchingCut.residualPeakSmallerThanStdNoise = [];
% 
% % Merge Templates after template matching
% P.mergeTemplates.merge = 1;
% P.mergeTemplates.upsampleFactor = 3;
% P.mergeTemplates.atCorrelation = .96;
% P.mergeTemplates.ifMaxDistSmaller = 2.5;      % in units std of noise
% P.mergeTemplates.ifMaxRelDistSmallerPercent = 25;
% 
% % Run BOTM Template Matching on whole data
% P.botm.run = 0;
% P.botm.Tf = 55;
% P.botm.cutLeft = 10;
% P.botm.prior = .0001;
