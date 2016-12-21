function [intervalIndices] = intervalAssign(spiketimes, intervals)
%%
%if size(intervals, 2) == 1
    tend = [intervals(:,1)' intervals(end,end)];
%else
%    tend = intervals(
%end

%ep = [1 4; 4 7; 7 12; 12 40]
%x = [0.2; 2.3; 11.9; 4.4; 10.3; 12; 10.4; 39.99; 40.0; 40.1]
%ep2 = [ep(:,1)' ep(end,end)];
spiketimes(spiketimes< tend(1) | spiketimes>tend(end)) = [];
%tend = ep2



findInterval = @(spiketime) max(1./(tend-spiketime));
[~, intervalIndices] = arrayfun(findInterval, spiketimes);
intervalIndices = intervalIndices - 1;
    %function idx = findInterval(spiketime)
    %    [~, idx] = min( tend - spiketime );
    %end

end