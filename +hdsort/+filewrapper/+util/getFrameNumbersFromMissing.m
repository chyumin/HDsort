function frameNo = getFrameNumbersFromMissing(firstFrame, lastFrame, mfv)

if isstruct(firstFrame)
    if isfield(firstFrame, 'missing_fns')
        frameNo = reconstrictFrames(firstFrame);
        return
    end
    
    
    S = firstFrame;
    firstFrame = S.first;
    lastFrame = S.last;
    mfv = [S.begin; S.length];
end

nGaps = size(mfv, 2);
n_missing = sum(mfv(2,:));

frameNo = firstFrame:(lastFrame-n_missing);
for g = 1:nGaps
    startIdx = find(frameNo == mfv(1,g));
    frameNo(startIdx:end) = frameNo(startIdx:end) + mfv(2,g);
end


    function numb =  reconstrictFrames(fn);
        numb = [];
        for ii = 1:numel(fn)
            f = fn(ii).first_fn:fn(ii).last_fn;
            
            if numel(fn(ii).missing_fns) > 1
                
                for jj = 1:size(fn(ii).missing_fns, 1)
                    begin_f = find(f  == fn(ii).missing_fns(jj, 1) );
                    end_f = begin_f + fn(ii).missing_fns(jj, 2);
                    
                    f(begin_f:end_f) = [];
                    
                    
                end
            end
            
            numb = [numb, f];
            
        end
        
    end

end