
function merged = merge(hdsort.epoch.) 
    merged = [];
    if isempty(hdsort.epoch.)
        return
    end
    
    hdsort.epoch. = sortrows(hdsort.epoch.);
    
    merged(1,:) = hdsort.epoch.(1,:);
    k = 1;
    for i=2:size(hdsort.epoch.,1)
        if hdsort.epoch.(i,1) <= merged(k,2)
            merged(k,:) = [min(hdsort.epoch.(i,1), merged(k,1)) max(hdsort.epoch.(i,2), merged(k,2))];
        else
            k = k+1;
            merged(k,:) = hdsort.epoch.(i,:);
        end
    end