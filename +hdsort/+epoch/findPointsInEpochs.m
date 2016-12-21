function fp = findPointsInEpochs(points, hdsort.epoch.)
    % finds all points in points that happen during an hdsort.epoch.in hdsort.epoch.
    pointPtr = 1;
    hdsort.epoch.tr = 1;
    
    points = sort(points);
    hdsort.epoch. = sortrows(hdsort.epoch.);
    
    nP = length(points);
    nE = size(hdsort.epoch.,1);
    
    fp = zeros(nP,1);
    while pointPtr<=nP && hdsort.epoch.tr<=nE
        if points(pointPtr) < hdsort.epoch.(hdsort.epoch.tr,1)
            % current point is before hdsort.epoch.
            pointPtr = pointPtr+1;
        elseif points(pointPtr) <= hdsort.epoch.(hdsort.epoch.tr,2)
            % current point is in current hdsort.epoch.
            fp(pointPtr) = hdsort.epoch.tr;
            pointPtr = pointPtr+1;
        else
            % current point is behind current hdsort.epoch.
            hdsort.epoch.tr = hdsort.epoch.tr+1;
        end        
    end        
end