function gdf = toGdfFrom2DCell(St, startTimes)

[nUnits, nFiles] = size(St);
if nargin < 2
    startTimes = zeros(1, nFiles);
end
assert(length(startTimes) == nFiles, 'Number of start times must be equal to number of spikes!')

gdf = [];
for f = 1:nFiles
    gdf_ = double(mysort.spiketrain.toGdf({St{:,f}}));
    l = size(gdf_, 1);
    gdf_ = gdf_ + [zeros(l,1), ones(l,1)*startTimes(f)];
    gdf = [gdf; gdf_];
end

end