function [T, tau, mMPmask] = tAlignOnUpsampleMean(T, varargin)
    P.upsample = 3;
    P.maxIter = 30;
    P.restrictToIdx = [];
    P.restrictToChannels = [];
    P.restrictToNMaximalValues = [];
    P.restrictMeanToItems = [];
    P.meanMaskMinLen = [];
    P.downsample = 1;
    P.maxIdx = [];
    P.initAlignment = '-';
    P.useMedian = 0;
    P.maxShiftPerIter = 3;
    P = mysort.hdsort.util.parseInputs(P, varargin, 'error');
    
    % Upsample all 
%     T = mysort.hdsort.util.resampleTensor(T, P.upsample, 1);
    if ~isempty(P.upsample) && P.upsample > 1
        disp('Upsampling Waveforms...')
        T = hdsort.waveforms.tResample(T, P.upsample, 1, 1);
    end
    if 0
        %%
        vT = hdsort.waveforms.t2v(T); 
        figure, hdsort.plot.vT(1:1000,:)')
    end
    % figure; hdsort.plot.squeeze(T))
    % Initial alignment on maxima
    if ~isempty(P.initAlignment)
        if isnumeric(P.initAlignment)
            tau_init = round(P.initAlignment*P.upsample);
            tau_init = tau_init(:);
            T = hdsort.waveforms.tShift(T, tau_init, 1);
        elseif strcmp(P.initAlignment, '+')
            [T tau_init] = hdsort.waveforms.tAlignOnMax(T, 'truncate', 1, ...
                'restrictToIdx', P.restrictToIdx, ...
                'restrictToChannels', P.restrictToChannels, 'maxIdx', P.maxIdx*P.upsample);
        elseif strcmp(P.initAlignment, '-')
            [T tau_init] = hdsort.waveforms.tAlignOnMax(-T, 'truncate', 1, ...
                'restrictToIdx', P.restrictToIdx, ...
                'restrictToChannels', P.restrictToChannels, 'maxIdx', P.maxIdx*P.upsample);
            T = -T;
        end
    else
        tau_init = zeros(size(T,3),1);
    end
    % figure; hdsort.plot.squeeze(T))
    % align on mean
    [T, tau_m, mMPmask] = hdsort.waveforms.tAlignOnMean(T, 'maxIter', P.maxIter, ...
        'useMedian', P.useMedian, ...
        'restrictToNMaximalValues', P.restrictToNMaximalValues*P.upsample,...
        'restrictToChannels', P.restrictToChannels,...
        'restrictMeanToItems', P.restrictMeanToItems, ...
        'meanMaskMinLen', P.meanMaskMinLen, ...
        'maxShiftPerIter', P.upsample*P.maxShiftPerIter);
    tau = tau_m+tau_init;
    
    if ~isempty(P.downsample) && P.downsample
        % Do this the brutal way:
        T = T(1:P.upsample:end,:,:);
        tau = tau/P.upsample;
        % Dont use resample again. This takes very long since it applies
        % anti aliasing. We do not need to do that, there shouldnt be any
        % high frequencies in the data, since we just upsampled it.
%         T = mysort.hdsort.util.resampleTensor(T, 1, P.upsample);
    end