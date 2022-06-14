function [RF, T] = receptiveFieldSVD(rf3d, varargin)
% Input:
%  - RF: receptive field as 2D-matrix.   
%  - varargin: [so far only 'debug']
%
% Output:
%  - out: structure containing all fitting parameters
%  - P: input preferences

P.debug = false;
%P.MaxIter = [];
P = hdsort.util.parseInputs(P, varargin);



rf3dDims = size(rf3d);
rf = reshape(rf3d, prod(rf3dDims(1:2)), rf3dDims(3));
[nPix, nTime] = size(rf);

[U,S,V] = svd(rf', 'econ');
% (time, space)
% U(:,1): time filter
% V(:,1): space filter

%
firstRFcomp = reshape(V, rf3dDims(1), rf3dDims(2), []);

RF = firstRFcomp(:,:,2);
T = U(:,1);

if P.debug
    implay(uint8(255*mat2gray(firstRFcomp)));%
    figure; hdsort.plot.U(:,1)
    figure; imagesc(firstRFcomp(:,:,2))
    figure; MG_Analysis.hdsort.plot.ectors(n)
end

end
