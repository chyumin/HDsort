function all_xc = computeXCorrs(X, cp, maxLag, hdsort.epoch., maxSamples)

    %
    % computes the pairwise xcorr functions of columns cp(:,1) and cp(:,2)
    % of matrix X in the hdsort.epoch. (parts of columns) given in hdsort.epoch..
    % if the length of all hdsort.epoch. is larger than maxSamples, only so many
    % hdsort.epoch. will be used, that at least maxSamples are contained.
    
    if isempty(cp)
        cp = mysort.hdsort.noise.computeChannelPairs4Channels(1:size(X,2));
    end
    
    nCP = size(cp,1);
%     if matlabpool('size')==0;
%         warning('No Matlabpool initialized. You can speed up computation by calling, e.g., matlabpool(n) with n being your computers core number');
%     end
    if nargin < 5
        maxSamples = size(X,1);
    end    
    if nargin < 4 || isempty(hdsort.epoch.)
        hdsort.epoch. = cell(1, nCP);
        for i=1:nCP
            hdsort.epoch.{i} = [1 min(size(X,1), maxSamples)];
        end
    elseif ~iscell(hdsort.epoch.)
        hdsort.epoch._ = hdsort.epoch.;
        hdsort.epoch. = cell(1, nCP);
        for i=1:nCP
            hdsort.epoch.{i} = hdsort.epoch._;
        end
    else
        assert(length(hdsort.epoch.) == nCP, 'There must be one hdsort.epoch.per channelPair!');
    end
    
    all_xc = zeros(2*maxLag+1, nCP);
    for i=1:nCP
        myEpochs = hdsort.epoch.{i};
        L = mysort.hdsort.epoch.length(myEpochs);
        totalEpochLength = sum(L);

        if totalEpochLength > maxSamples
            % randomy permute hdsort.epoch.
            rp = randperm(size(myEpochs,1));
            cs = cumsum(L(rp));
            lastIdx = find(cs>maxSamples,1);
            myEpochs = myEpochs(rp(1:lastIdx),:);
    %         L = L(rp(1:lastIdx));
            totalEpochLength = cs(lastIdx);
        end        
        for k=1:size(myEpochs,1)
            % get all data for this hdsort.epoch.
            x = X(myEpochs(k,1):myEpochs(k,2), cp(i,:));
            % compute channel pairs
            all_xc(:,i) = all_xc(:,i)+xcorr(x(:,1), x(:,2), maxLag, 'none');        
        end

        % take only those who are needed and normalize
        all_xc(:,i) = all_xc(:,i)/totalEpochLength;        
    end
    
    
    
    
    
%% GRAVEYARD
% old way of computing xcorrs:

    %%% WARNING! %%%
    % This implementation assumes that it is faster to compute all cross
    % correlations between the whole set of channels either in ep(:,1) OR
    % in ep(:,2) rather than computing only those pairs that are actually
    % needed (by some for loop).
    %%%
%     % replace the channel index into X by the relative channel idx
%     uChans = unique(cp(:));
%     [a cpidx] = ismember(cp, uChans);
% 
%     xc = zeros(2*maxLag+1, length(uChans)^2);
%     for k=1:size(hdsort.epoch.,1)
%         % get all data for this hdsort.epoch.
%         x = X(hdsort.epoch.(k,1):hdsort.epoch.(k,2), uChans);
%         % compute channel pairs
%         xc = xc+xcorr(x, maxLag, 'none');        
%     end
%     nCp = size(cp,1);
%     % compute the index into xc for each of the cp
%     xcidx = (cpidx(:,1)-1)*length(uChans) + cpidx(:,2);
%     
%     % take only those who are needed and normalize
%     xc = xc(:, xcidx)/totalEpochLength;
