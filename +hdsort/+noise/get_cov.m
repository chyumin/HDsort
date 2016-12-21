

%% Compute hdsort.noise.hdsort.epoch. common with n nearest electrodes
NEN = NE;
for ch = 1:nC
    % find k nearest electrodes
    neighbors = mysort.mea.nearestElectrodes(channel_x, channel_y, ch, 6);
    for n = 1:length(neighbors)
        NEN{ch} = mysort.hdsort.epoch.intersect(NEN{ch}, NE{neighbors(n)});
    end
end


%%
nXc = (nC*(nC-1))/2;
D = zeros(1,nXc);count = 1;
Tf = 30; maxLag = Tf-1;
autocovs = zeros(nC, 2*Tf-1);
autocov_norms = zeros(nC,1);
xcovs = zeros(nXc, 2*Tf-1);
IDX = zeros(nXc,2);


for ch1 = 1:nC
    p1 = [channel_x(ch1) channel_y(ch1)];
    [ac lc] = mysort.hdsort.util.xcorr_in_hdsort.epoch.(Y(ch1,:), NE{ch1}, maxLag, maxLag);
    autocovs(ch1,:) = ac;
    autocov_norms(ch1) = lc;
    
    for ch2 = ch1+1:nC 
        p2 = [channel_x(ch2) channel_y(ch2)];
        d = norm(p1-p2);
        IDX(count,:) = [ch1 ch2];
        D(count) = d; 
        if d < 40
            cNE = mysort.hdsort.epoch.intersect(NE{ch1},NE{ch2});
            cNE = mysort.hdsort.epoch.removeShort(cNE, 10*Tf);
            cNE = cNE(1:3:end,:);
            [xc lc] = mysort.hdsort.util.xcorr_in_hdsort.epoch.(Y([ch1 ch2],:), cNE, maxLag, maxLag);
        else
            xc = zeros(1, 2*Tf-1);
        end
        xcovs(count,:) = xc;   
        count=count+1;
        disp(count)
    end
end
% %%
% save('C:\LocalData\Michele\marching_square_buffer_covs', 'D', 'Tf', 'autocovs', ...
%      'autocov_norms', 'D', 'xcovs', 'IDX', 'smad', 'NE');
% %%
% hdsort.noise.dsort.epoch1 = mysort.hdsort.epoch.flip(spikehdsort.epoch.1, size(X,2));
% [smad2, spikehdsort.epoch.2] = ana.douglas.estimateSigma(X(k2,:), Tf, thr1);
% hdsort.noise.dsort.epoch2 = mysort.hdsort.epoch.flip(spikehdsort.epoch.2, size(X,2));
% 
% spikehdsort.epoch.1_idx = mysort.hdsort.epoch.toIdx(spikehdsort.epoch.1);
% spikehdsort.epoch.2_idx = mysort.hdsort.epoch.toIdx(spikehdsort.epoch.2);
% hdsort.noise.dsort.epoch1_idx = mysort.hdsort.epoch.toIdx(hdsort.noise.dsort.epoch1);
% hdsort.noise.dsort.epoch2_idx = mysort.hdsort.epoch.toIdx(hdsort.noise.dsort.epoch2);
% commonhdsort.noise.dsort.epoch = mysort.hdsort.epoch.intersect(hdsort.noise.dsort.epoch1,hdsort.noise.dsort.epoch2);
% commonhdsort.noise.dsort.epoch_idx = mysort.hdsort.epoch.toIdx(commonhdsort.noise.dsort.epoch);
% unionspikehdsort.epoch.  = mysort.hdsort.epoch.merge( [spikehdsort.epoch.1; spikehdsort.epoch.2]);
% unionspikehdsort.epoch._idx = mysort.hdsort.epoch.toIdx(unionspikehdsort.epoch.);