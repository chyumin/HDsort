
function R = compare(st1, st2, maxJitter, maxShift, useparfor)
if nargin < 5
    useparfor = true;
end

n1 = length(st1);
n2 = length(st2);
R.N1 = zeros(n1,1);
R.N2 = zeros(n2,1);
R.TP = zeros(n1, n2);
R.FP = zeros(n1, n2);
R.FN = zeros(n1, n2);
R.ALI = {};
R.CL1 = {};
R.CL2 = {};
R.readme = 'Contains the result data of a comparison between two sets of spike trains. It is assumed, that the first set is some kind of ground truth information, meaning the true spike trains. The second set might be the result of a spike sorting. Then, for every combination between the a true spike train (rows of the matrices) and a sorted one (columns), a comparison is made.';
if isempty(maxShift)
    maxShift = 0;
end
shiftRange = -maxShift:floor(maxJitter/5):maxShift;
nullidx = find(shiftRange==0);
% Make sure zero shift is included and at the beginning of the shift
% range. this will ensure that max operations on the tp array later
% will favour the 0 shift if more than this have the same tp number
if ~isempty(nullidx)
    shiftRange(2:nullidx) = shiftRange(1:nullidx-1);
    shiftRange(1) = 0;
else
    shiftRange = [0 shiftRange];
end

for i=1:n1
    disp([num2str(i) ' of ' num2str(n1) ' computed...'])
    
    R.N1(i) = length(st1{i});
    
    if ~useparfor
        for j=1:n2
            R.N2(j) = length(st2{j});
            for shiftIdx = 1:length(shiftRange)
                shift = shiftRange(shiftIdx);
                [tp(shiftIdx), fp(shiftIdx), fn(shiftIdx), alignment{shiftIdx}, cla1{shiftIdx}, cla2{shiftIdx}] = ...
                    hdsort.spiketrain.compareTwo(st1{i}, st2{j}+shift, maxJitter);
            end
            [mTp mTpIdx] = max(tp);
            R.shift(i,j) = shiftRange(mTpIdx);
            R.TP(i,j) = tp(mTpIdx);
            R.FP(i,j) = fp(mTpIdx);
            R.FN(i,j) = fn(mTpIdx);
            R.ALI{i,j} = alignment{mTpIdx};
            R.CL1{i,j} = cla1{mTpIdx};
            R.CL2{i,j} = cla2{mTpIdx};
        end
        
    else
        
        R_ = {};
        parfor j=1:n2
            R_{j}.N2 = length(st2{j});
            
            L = length(shiftRange);
            tp = ones(1,L); fp = ones(1,L); fn = ones(1,L); alignment = {}; cla1 = {}; cla2 = {};
            for shiftIdx = 1:length(shiftRange)
                shift = shiftRange(shiftIdx);
                
                [tp(shiftIdx), fp(shiftIdx), fn(shiftIdx), alignment{shiftIdx}, cla1{shiftIdx}, cla2{shiftIdx}] = ...
                    hdsort.spiketrain.compareTwo(st1{i}, st2{j}+shift, maxJitter);
            end
            [mTp mTpIdx] = max(tp);
            R_{j}.shift = shiftRange(mTpIdx);
            R_{j}.TP = tp(mTpIdx);
            R_{j}.FP = fp(mTpIdx);
            R_{j}.FN = fn(mTpIdx);
            R_{j}.ALI = alignment{mTpIdx};
            R_{j}.CL1 = cla1{mTpIdx};
            R_{j}.CL2 = cla2{mTpIdx};
        end
        
        for j=1:n2
            R.N2(j) = R_{j}.N2;
            R.shift(i,j) = R_{j}.shift;
            R.TP(i,j) = R_{j}.TP;
            R.FP(i,j) = R_{j}.FP;
            R.FN(i,j) = R_{j}.FN;
            R.ALI{i,j} = R_{j}.ALI;
            R.CL1{i,j} = R_{j}.CL1;
            R.CL2{i,j} = R_{j}.CL2;
        end
    end
    
end
end