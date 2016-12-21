function P = boxhdsort.plot.data, groupIdx, varargin)
% Data a NxS matrix containing N datapoints and S dimensions (subgroups)
% each.
P.fh = [];
P.ah = [];
P.figureName = 'Boxhdsort.plot.;
P.ylabel = ''
P.displayN = false;

P.subgroupnames = {};
P.groupnames = {};
P = mysort.hdsort.util.parseInputs(P, varargin, 'error');

if size(data, 1) == 1
    data = data';
end
    
if isempty(P.fh) & isempty(P.ah)
    P.fh = mysort.hdsort.plot.figure('name',P.figureName);
elseif ~isempty(P.fh) & isempty(P.ah)
    figure(P.fh);
else
    P.ah = axes();
end

groups = unique(groupIdx);
nGroups = length(groups);
[nDataPoints nSubGroups] = size(data);
%nSubGroups = nSubGroups_/nGroups;

%% Create labels for groups:
if isempty(P.groupnames)
    for g = 1:nGroups
        P.groupnames{g} = num2str(groups(g));
    end
end
grouplabels = repmat(P.groupnames,1,nSubGroups);

%% Create labels for subgroups:
if isempty(P.subgroupnames)
    for sg = 1:nSubGroups
        P.subgroupnames{sg} = num2str(sg);
    end
end

subgrouplabels = [];
for sg = 1:nSubGroups
    subgrouplabels = [subgrouplabels, repmat({P.subgroupnames{sg}},1,nGroups)];
end

col = mysort.hdsort.plot.vectorColor(1:nSubGroups);

%% Reshape data:
data_ = zeros(nDataPoints, nGroups*nSubGroups)*NaN;
for g = 1:nGroups
    dataGroup = data(groupIdx == groups(g), :);
    [L nSubGroups_] = size(dataGroup);
    assert(nSubGroups==nSubGroups_, 'error')
    
    for sg = 1:nSubGroups
        ii =(sg-1)*nGroups + g;
        data_(1:L,ii) = dataGroup(1:L,sg);
    end
end

if nSubGroups == 1
    boxhdsort.plot.data_, grouplabels);
else
    P.boxhdsort.plot.= boxhdsort.plot.data_,{grouplabels,subgrouplabels}, 'colors', col,'factorgap',[5 2],'labelverbosity','minor');
end
ylabel(P.ylabel);

if P.displayN
    M = nanmean(data_);
    N = sum(~isnan(data_));
    y = max(data_(:));
    for k1 = 1:size(data_,2)
        text(k1, y, sprintf('N = %d', N(k1)), 'FontSize',8);
    end
end

end
