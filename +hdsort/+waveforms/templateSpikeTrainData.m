function [xticks Y] = templateSpikeTrainData(T, gdf, cutLeft)
    Tf = size(T,1);
    nC = size(T,2);
    nT = size(T,3);
    
    if nargin >4
        gdf = gdf(gdf(:,2)>=t1 & gdf(:,2)<=t2,:);
    end
    if isempty(gdf)
        xticks = [];
        Y = [];
        return
    end
    
    % add amplitudes if necessary
    if size(gdf,2) == 2
        gdf = [gdf ones(size(gdf,1),1)];
    end
    
    hdsort.epoch. = [gdf(:,2)-cutLeft gdf(:,2)-cutLeft + Tf-1];
    hdsort.epoch.M = mysort.hdsort.epoch.merge(hdsort.epoch.);
    hdsort.epoch.MLength = mysort.hdsort.epoch.length(hdsort.epoch.M);
    
    dataL = sum(hdsort.epoch.MLength) + size(hdsort.epoch.M,1);
    hdsort.epoch.MstartInData = cumsum(hdsort.epoch.MLength+1)';
    hdsort.epoch.MstartInData = [0 hdsort.epoch.MstartInData(1:end-1)]+1;
    
    xticks = zeros(1, dataL);
    
    for i=1:size(hdsort.epoch.M,1)
        s1 = hdsort.epoch.MstartInData(i);
        s2 = s1+hdsort.epoch.MLength(i)-1;
        xticks(s1:s2+1) = [hdsort.epoch.M(i,1):hdsort.epoch.M(i,2) nan];
    end
    
    Y = zeros(nC, dataL);
    
    myEpochIdx = 1;
    for i=1:size(gdf,1)
        if gdf(i,2) > hdsort.epoch.M(myEpochIdx,2)
            myEpochIdx = myEpochIdx +1;
        end
        
        offsetInEpoch = gdf(i,2)-cutLeft - hdsort.epoch.M(myEpochIdx,1);
        
        s1 = hdsort.epoch.MstartInData(myEpochIdx)+offsetInEpoch;
        s2 = s1+Tf-1;
        if s1>0 && s2 <= size(Y,2)
            Y(:, s1:s2) = Y(:, s1:s2)+gdf(i,3)*[squeeze(T(:,:,gdf(i,1)))'];
        end
%             warning()
    end
    
    
    
    

