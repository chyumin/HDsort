function [spike_trains uIDs] = gdf2cell(tGDF)
spike_trains = {};
if isempty(tGDF)
    return
end
%tIDs = unique(tGDF(:,2));
uIDs  = unique(tGDF(:,1));
%nT = length(tIDs);
nU = length(uIDs);
spike_trains = cell(1, nU);

for unit = 1:nU
    spike_trains{unit} = tGDF(tGDF(:,1) == uIDs(unit), 2);
end
