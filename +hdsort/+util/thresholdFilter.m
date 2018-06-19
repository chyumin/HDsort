function [SELECTION_INDEX, or_idx, xor_idx] = thresholdFilter(thresholds, S, printFlag) 
    if nargin < 3
        printFlag = 0;
    end
    % Selects elements of structure array S given thresholds for its member
    % fields in tresholds
    nThr = size(thresholds,1);
    nUTotal = length(S);
    indices = zeros(nUTotal, nThr);
    for ti = 1:nThr
        if isa(thresholds{ti, 2}, 'function_handle')
            indices(:, ti) = thresholds{ti, 2}([S.(thresholds{ti,1})]);
        elseif strcmp(thresholds{ti, 2}, '>')
            indices(:, ti) = [S.(thresholds{ti,1})] > thresholds{ti,3};
        elseif strcmp(thresholds{ti, 2}, '<')
            indices(:, ti) = [S.(thresholds{ti,1})] < thresholds{ti,3};
        elseif strcmp(thresholds{ti, 2}, '==')
            indices(:, ti) = [S.(thresholds{ti,1})] == thresholds{ti,3};
        elseif strcmp(thresholds{ti, 2}, '~=')
            indices(:, ti) = [S.(thresholds{ti,1})] ~= thresholds{ti,3};            
        elseif strcmp(thresholds{ti, 2}, '~') && isempty(thresholds{ti,3})
            indices(:, ti) = ~isempty([S.(thresholds{ti,1})]);  
        elseif strcmp(thresholds{ti, 2}, 'is') && isempty(thresholds{ti,3})
            indices(:, ti) = isempty([S.(thresholds{ti,1})]);              
        end
    end

    %
    SELECTION_INDEX = all(indices, 2);
    SELECTION_SET = find(SELECTION_INDEX);

    or_idx  = sum(indices==0,1);
    xor_idx = sum(indices==0 & repmat(sum(indices==0,2) == 1, [1, nThr]));

    if printFlag
        if isempty(thresholds)
            fprintf('##### %d (ALL) Neurons preselected - No thresholds provided\n', length(SELECTION_SET));
        else
            fprintf('##### %d Neurons preselected\n', length(SELECTION_SET));
        end
    end