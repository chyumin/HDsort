
function [ST, classes, ori_classes, AMPS] = fromGdf(gdf)
ST = {}; classes = [];
if isempty(gdf)
    return
end
assert(size(gdf,2)>=2, 'A gdf has to be a matrix with two or more colums!');


ori_classes = unique(gdf(:,1));
assert(length(ori_classes) < 3000, 'This function should not be used for more than 3000 neurons!');
classes = ori_classes;
zer = find(classes==0,1);
if ~isempty(zer)
    gdf(:,1) = gdf(:,1) +1;
    classes = classes+1;
end
ST = cell(1, length(classes));
myIdx = [];

if size(gdf,2) >= 4
    AMPS = cell(1, length(classes));
    for ii = 1:length(classes)
        myClass = classes(ii);
        myIdx = gdf(:,1)==myClass;
        
        x = gdf(myIdx,[2,4]);
        ST{ii} = x(:,1);
        AMPS{ii} = x(:,2);
    end
else
    AMPS = [];
    for ii = 1:length(classes)
        myClass = classes(ii);
        myIdx = gdf(:,1)==myClass;
        
        ST{ii} = gdf(myIdx,2);
    end
end

end