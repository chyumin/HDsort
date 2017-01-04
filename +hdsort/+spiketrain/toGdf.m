
function gdf = toGdf(St, IDs)
    gdf = [];
    if isempty(St)
        return
    end
    if nargin>1
        assert(length(St) == length(IDs), 'There must be one ID per spike train!');
    else
        IDs = 1:length(St);
    end
%    assert(size(gdf,2)==2, 'A gdf has to be a matrix with exactly two colums!');
    for i=1:length(St)
        if ~isempty(St{i})
            gdf = [gdf; [ones(length(St{i}),1)*IDs(i) hdsort.util.toColVec(St{i})]];
        end
    end
    if isempty(gdf)
        return
    end
    gdf = sortrows(gdf,2);
end