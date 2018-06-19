function frameNo = getFrameNumbersFromMissing(missingFrameNumbers)

if numel(missingFrameNumbers) > 1
    frameNo = [];
    for ii = 1:numel(missingFrameNumbers)
        frameNo = [frameNo, hdsort.file.util.getFrameNumbersFromMissing(missingFrameNumbers(ii))];
    end
    return;
end

frameNo = missingFrameNumbers.first:missingFrameNumbers.last;
idx = [];
for ii = 1:numel(missingFrameNumbers.begin)
    i0 = find(frameNo == missingFrameNumbers.begin(ii));
    idx = [idx, i0:i0+missingFrameNumbers.length(ii)-1];
end
frameNo(idx) = [];

end
