function [FR, x] = toFiringRate(ST, window, sampling, timePeriod)
    % converts spike trains to vectors of firing rates
    % Input:
    %  ST       - cell array of spike Trains in SECONDS !
    %  window   - width of kernel for smoothing in seconds
    %  sampling - sampling (delta t) of the returned firing rate in
    %             seconds!
    %  timePeriod  - vector with two elements giving the temporal length of
    %  (optional)    the resulting firing rate vectors as start and stop
    %                time (in seconds)
    %
    % Output:
    %  FR       - time series of firing rates if ST is a vector
    %             cell array of time series if ST is a cell array
    
    x = [];
%     if isempty(ST)
%         FR = ST;
%         return
%     end

    
    if iscell(ST)
        if nargin < 4
            allst = cat(1, ST{:});
            allst = allst(:);
            timePeriod = [min(allst) max(allst)];
        end
        FR = cell(size(ST));
        for i = 1:size(ST,1)
            for j = 1:size(ST,2)
                [FR{i,j}, x] = mysort.spiketrain.toFiringRate(ST{i,j}, window, sampling, timePeriod);
            end
        end
        return
    else
        if nargin < 4
            timePeriod = [min(ST(:)) max(ST(:))];
        end        
        timePeriod = timePeriod(:);
        assert(length(timePeriod) == 2, 'must be 2 element vector');
        minTime = timePeriod(1);
        maxTime = timePeriod(2);
        nBins = ceil((maxTime-minTime)/sampling) + 1;
        x = (0:nBins-1)*sampling + minTime;
        if isempty(ST)
            FR = zeros(length(x), 1);
        else
            FR = histc(ST, x);
            if ~isempty(window)
                windowInSamples = ceil(window/sampling);
                tauRange = -2*windowInSamples:2*windowInSamples;
                kernel = normpdf(tauRange, 0, windowInSamples/2);
                
                % This normalizes the kernel to an integral of 1 in respect
                % to the temporal sampling 
                kernelI = trapz(tauRange*sampling, kernel); % normalization?
                kernel = kernel/kernelI;
%                 kernelI2 = trapz(tauRange*sampling, kernel);

                FR = conv(FR, kernel, 'same');
            end
        end
        FR = FR(:);
    end
    
    
    
    
    
    