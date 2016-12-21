function [centers, N_centers] = createHexagonalGrid(N_circles, varargin)

P.debug = false;
P.distances = 1.0;
P.center = [0,0];

P = hdsort.util.parseInputs(P, varargin, 'error');




radius = N_circles + 0.1;
Rad3Over2 = sqrt(3) / 2;
[X Y] = meshgrid(1:1:11);
n = size(X,1);
X = Rad3Over2 * X;
Y = Y + repmat(0.5*(1:n), size(Y,1), 1);
centers = [X(:), Y(:)];
mid_center = median(centers);
centers(pdist2(centers, mid_center) > radius, :) = [];

shift = mid_center - P.center/P.distances;

centers = centers - repmat(shift, size(centers,1), 1);
centers = centers * P.distances;

if P.debug
    figure; scatter(centers(:, 1),centers(:,2))
end

N_centers = size(centers, 1);

end