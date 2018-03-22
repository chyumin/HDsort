function [filteredFP, P] = filterFootprints(FP, varargin)

%% Design waveform FIR filter
%P.subtractMeanOverAllChannels = true; % already done before in zero offset function!
P.hpf = 300;
P.lpf = 6000;
P.fir_filterOrder = 110;
P.frameRate = 20000;
P.doZeroPad = false;
P = util.parseInputs(P, varargin, 'error');

%% Create Filter:
P.fir_filter = mysort.mea.filter_design_fir(P.hpf, P.lpf, P.frameRate, P.fir_filterOrder);

%%
[nTf, nCh, nU] = size(FP);

if P.doZeroPad
    
    P.zeroPad = ceil(numel(P.fir_filter)*1.5);
    z = ones(P.zeroPad, 1);
    
    filteredFP = zeros(nTf, nCh, nU);
    for ch = 1:nCh
        for uu = 1:nU
            
            f_ = conv([z*FP(1,ch,uu); FP(:,ch,uu); z*FP(end,ch,uu)], P.fir_filter, 'same');
            filteredFP(:, ch, uu) = f_(P.zeroPad+1:end-P.zeroPad);
        end
    end
    
else
    
    L = P.fir_filterOrder;
    filteredFP = zeros(nTf-2*L, nCh, nU);
    
    for ch = 1:nCh
        for uu = 1:nU
            %filteredFP(:, ch, uu)  = conv(FP(:,ch,uu), P.fir_filter, 'same');
            f = conv(FP(:,ch,uu), P.fir_filter, 'same');
            filteredFP(:, ch, uu) = f(1+L:end-L);
            
        end
    end
    
end

