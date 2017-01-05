function [R, P] = processLocalSortings(dpath, runName, elGroupNumbers, elGroupIndices, varargin)
P.groupPaths = [];
P.newPostProcFunc = '';
P = hdsort.util.parseInputs(P, varargin, 'error');

%% Load GDFs
nGroups = length(elGroupNumbers);
for ii=1:nGroups
    group_paths       = fullfile(dpath, 'groups', sprintf('group%04d', ii));
    
    if iscell(P.groupPaths)
        assert(length(P.groupPaths) == nGroups, 'Something awefully horrible has happended!');
        subdpath{ii} = fullfile(dpath, P.groupPaths{ii});
    else
        assert(exist(group_paths, 'dir')==7, 'Something awefully horrible has happended!');
        subdpath{ii} = group_paths;
    end 
end

G = struct('sortedElectrodes', {}, 'gdf', {}, 'units', {});
nUnitsPerLocalSorting = zeros(1, nGroups);

parfor ii = 1:nGroups
    spdet = load(fullfile(subdpath{ii}, [runName '.030spikes_det_merged.mat']));
    clust = load(fullfile(subdpath{ii}, [runName '.110clusters_meanshift_merged.mat']));
    botmm = load(fullfile(subdpath{ii}, [runName '.100botm_matching.mat']));
    
    % This line will not work if P.spikeCutting.maxSpikes was set to
    % low in the call for the spike sorting:
    assert(length(clust.clusteringMerged.ids) == length(spdet.spikeDetectionMerged.ts), ...
        ['Not all spikes that were detected were cut (and thus also not matched). This is' ...
        'necessary for the postprocessing. Set P.spikeCutting.maxSpikes higher if you want' ...
        'all spikes to be cut-and-matched.']);
    
    if ~isfield(botmm.clusteringMatched, 'amps')
        botmm.clusteringMatched.amps = 0*clust.clusteringMerged.ids;
        botmm.clusteringMatched.ampsChIdx = 0*clust.clusteringMerged.ids + 1;
    end
    gdf = [clust.clusteringMerged.ids spdet.spikeDetectionMerged.ts ...
        botmm.clusteringMatched.amps elGroupIndices{ii}(botmm.clusteringMatched.ampsChIdx)];
    units = unique(gdf(:,1));
    
    G(ii).sortedElectrodes = elGroupIndices{ii};
    G(ii).gdf = gdf;
    G(ii).units = units;
    nUnitsPerLocalSorting(ii) = length(units);
end


fprintf('Found Units in local sortings: %d\n', sum(nUnitsPerLocalSorting));

%% Load Templates
nT = 0;
for ii=1:nGroups
    fnames = dir(fullfile(subdpath{ii}, [runName '_templates.mat']));
    G(ii).templates = load(fullfile(subdpath{ii}, fnames.name));
    nT = nT + size(G(ii).templates.wfs,3);
end
fprintf('Found Templates local sortings: %d\n', nT);

%% Load Noise
Ns = cell(nGroups,1);
for ii=1:nGroups
    S = load(fullfile(subdpath{ii}, [runName '.060cov.mat']));
    Ns{ii} = S.noise.CestS.CCol;
end
meanNoiseStd = mean(sqrt(diag(Ns{1})));

%% Merge double Templates
for ii=1:length(G)
    nUnitsInGDF      = length(unique(G(ii).gdf(:,1)));
    
    if ~isempty(G(ii).templates.wfs)
        nUnitsInTempates = size(G(ii).templates.wfs,3);
    else
        nUnitsInTempates = 0;
    end
    str_ = sprintf('ii: %d, nGdf: %d, nWfs: %d, dpath: %s, rn: %s', ii, nUnitsInGDF, nUnitsInTempates, dpath, runName);
    assert(nUnitsInGDF == nUnitsInTempates, ['must be identical - ' str_])
end

fprintf('Merging double Templates...');
if strcmp(P.newPostProcFunc, '')
    G = hdsort.leg.mergeLocalSortings(G, meanNoiseStd);
else
    disp(['Postprocessing function: ' P.newPostProcFunc]);
    fh = str2func(P.newPostProcFunc);
    G = fh(G, meanNoiseStd);
end

%     G = hdsort.leg.combineLocalSortings(G, meanNoiseStd);
for ii=1:length(G)
    if ~isempty(G(ii).gdf(:,1))
        assert(length(unique(G(ii).gdf(:,1))) == size(G(ii).templates.wfs,3), ['must be identical - ' str_])
    else
        assert(isempty(G(ii).templates.wfs), ['must be identical - ' str_])
    end
end
g_file = fullfile(dpath, 'G_struct.mat');
save(g_file, 'G', 'meanNoiseStd', '-v7.3');

%% Extract final gdf and templates
disp('Extracting final gdf and templates...');
[R, P_] = hdsort.leg.computeFinalGDFFromMergedLocalSortings(G);

fprintf('Templates after merging: %d\n', size(R.T_merged,3) );
R.gdf_merged = sortrows(R.gdf_merged,2);

if ~isempty(R.gdf_discarded)
    R.gdf_discarded = sortrows(R.gdf_discarded, 2);
end
R.G = G;

