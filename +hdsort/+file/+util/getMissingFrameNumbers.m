function missingFrames = getMissingFrameNumbers(frame_numbers)
assert(size(frame_numbers, 1) == 1, 'Input must be a row vector!')

df = [diff(frame_numbers)-1 0];
assert(~any(df<0), 'Frame numbers are not monotonically increasing!')

idx = find(df);

missingFrames.begin = frame_numbers(idx)+1;
%missingFrames.end = frame_numbers(idx)+1+df(idx);
missingFrames.length = df(idx);
missingFrames.n = numel(missingFrames.begin);

missingFrames.first = frame_numbers(1);
missingFrames.last  = frame_numbers(end);