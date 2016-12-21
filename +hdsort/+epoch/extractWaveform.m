
function spikesX = extractWaveform(X, hdsort.epoch.)  
    spikesX =  [];
    if isempty(hdsort.epoch.)
        return
    end
    assert(size(hdsort.epoch.,2) == 2, 'hdsort.epoch. must be a two column matrix!');
    nSpikes = size(hdsort.epoch.,1);
    nC = size(X,1);
    Tf = hdsort.epoch.(1,2)-hdsort.epoch.(1,1)+1;
    spikesX = zeros(nSpikes, nC*Tf);

    
    for i=1:size(hdsort.epoch.,1)
        startSample = hdsort.epoch.(i,1);
        endSample = hdsort.epoch.(i,2);       
        zL = 0;
        if startSample<=0 
            zL = -startSample+1;
            startSample=1;
        end
        zR = 0;
        if endSample>size(X,2)
            zR = endSample-size(X,2);
            endSample = size(X,2);
        end
        for c=1:nC
            spikesX(i,(1:Tf) + (c-1)*Tf) = ...
                [zeros(1,zL) X(c, startSample:endSample) zeros(1,zR)];
        end
    end
          