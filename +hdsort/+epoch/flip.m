

function flipped = flip(hdsort.epoch.,len) 
% enough comment!!!

    flipped = [];
    if isempty(hdsort.epoch.)
        if nargin>1
            flipped = [1 len];
        end
        return        
    end

    if hdsort.epoch.(1,1) > 1
        flipped(1,:) = [1 hdsort.epoch.(1,1)-1];
    end
    flipped = [ flipped;
            [hdsort.epoch.(1:end-1,2)+1 hdsort.epoch.(2:end,1)-1] ];
    
    % Remove hdsort.epoch. of negative of zero length
    flipped(flipped(:,2)-flipped(:,1)<1,:) = [];
        
    if (nargin>1) && hdsort.epoch.(end,2)<len
        flipped = [flipped; 
            [hdsort.epoch.(end,2)+1 len] ];
    end    
    
