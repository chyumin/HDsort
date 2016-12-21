function [Ts tau] = vComputeSubsampleShiftedVersions(vT, nC, upsample)
    % Computes the upsampled version of T shifts them (up)sample wise.
    % Then, each possible subsample shifted version is downsampled again
    % and returned in T. 
    % Upsampling is done with Matlabs resample function, but also
    % hdsort.waveforms.mSincfun would be possible
    %
    % Input:
    %    T        - 
    %    nC       -
    %    upsample - 
    %
    % Output:
    %    Ts - 4 dimensional tensor, time x channels x templates x subshifts
    %  tau  - 
    tT = hdsort.waveforms.v2t(vT, nC);
    [tT tau] = hdsort.waveforms.tComputeSubsampleShiftedVersions(tT, upsample);
    Ts = hdsort.waveforms.t2v(tT);
    
    