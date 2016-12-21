function P = gscatter(X, Y, labels, varargin)
%gscatter(AMP, Gt, alignmentLabel, col)

% Data a NxS matrix containing N datapoints and S dimensions (subgroups)
% each.
P.fh = [];
P.ah = [];
P.ylabel = '';
P.xlabel = '';
P.title = '';
P.markerSize = [];
P.color = [];
P = hdsort.util.parseInputs(P, varargin, 'error');


uLabels = unique(labels);
nGroups = numel(uLabels);

if isempty(P.fh) & isempty(P.ah)
    P.fh = hdsort.plot.figure('name',P.figureName);
    P.ah = axes();
elseif ~isempty(P.fh) & isempty(P.ah)
    figure(P.fh);
    P.ah = axes();
else
    axes(P.ah);
end
hold on;

if isempty(P.color)
    P.color = hdsort.plot.vectorColor(1:nGroups);
end

for g = 1:nGroups
    idx = labels == uLabels(g);
    scatter(P.ah, X(idx), Y(idx), P.markerSize, P.color(g,:) );
end

title(P.title);
xlabel(P.xlabel);
ylabel(P.ylabel);