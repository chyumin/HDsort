function [out P] = mexicanHatFit2(RF, varargin)
% Input:
%  - RF: receptive field as 2D-matrix.   
%  - varargin: [so far only 'debug']
%
% Output:
%  - out: structure containing all fitting parameters
%  - P: input preferences

P.debug = false;
P.MaxIter = [];
P = hdsort.util.parseInputs(P, varargin);

xrange = 1:size(RF,1);
yrange = 1:size(RF,2);

out.c = median(RF(:));
RF = RF - out.c;

if P.debug
    figure; imagesc(RF);
    figure;
end

EH = @(params) E(params, RF, xrange, yrange);

[max_x max_y] = hdsort.util.matrixArgMax(abs(RF));
ampmax = RF(max_x,max_y);
%params0 = [max_x max_y 10*rand(1) 10*rand(1) 10*rand(1) ampmax];
params0 = [max_x max_y var(sum(RF)) var(sum(RF')) rand(1)*2*pi ampmax];


if ~isempty(P.MaxIter)
    [paramsopt, fopt, exitflag] = fminsearch(EH, params0, optimset('MaxIter', P.MaxIter));
else
    [paramsopt, fopt, exitflag] = fminsearch(EH, params0); 
end


    function [mu C Amplitude U l ] = interpretParameters(params)
        mu = params(1:2);
        
        %% Constrain the difference between l_pos(1) and l_pos(2) to a ratio of epsilon:
            epsilon = 0.1;
            lower_bound = log(epsilon)/log(10);
            upper_bound = abs(lower_bound);
            
            l(1) = abs(params(3));
            l(2) = abs(params(3));
            %l(2) =  10^( (upper_bound-lower_bound)/( 1+exp(- params(4) ) ) + lower_bound) * l(1) ;
        
        %l_neg = abs(params(5)) * l_pos;
        
        u11 = sin(params(5));
        u22 = sqrt(1 - u11^2);
        u1 = [u11 u22];
        u2 = [-u1(2) u1(1)];
        
        U = [u1; u2];
        L = [l(1) 0; 0 l(2)];
        C = U*L*U';
       
        assert( ~any([C(1,1), C(2,2)] < 0), 'Cov matrices not well defined!')
        
        Amplitude = params(6);
    end

    function e = E(params, RF, xrange, yrange)
        
        [mu C Amplitude] = interpretParameters(params);
        
        RF_fit = RFmodel(xrange, yrange, mu, C, Amplitude);
        
        e = norm(RF_fit(:) - RF(:));
        
        if P.debug
            Amplitude
            %e
            subplot.2,1,1)
            %imagesc(reshape(RF_fit(:), size(RF,1), size(RF,2)) )
            imagesc(RF_fit)
            subplot.2,1,2)
            %imagesc(reshape(RF(:), size(RF,1), size(RF,2)) )
            imagesc(RF)
            pause(0.1)
        end
    end

    function rf = RFmodel(xrange, yrange, mu, Cov, Amplitude)
        [X Y] = meshgrid(xrange, yrange);
        XX = [X(:) Y(:)];
        rf = mvnpdf(XX, mu, Cov)*Amplitude;
        rf = reshape(rf, length(yrange), length(xrange))';
    end

[mu_opt C_opt a_opt U_opt l_opt] = interpretParameters(paramsopt);
out.RF = RFmodel(xrange, yrange, mu_opt, C_opt, a_opt);

if P.debug
    figure;
    subplot.2,1,1); imagesc(RF);
    subplot.2,1,2); imagesc(out.RF)
end

out.mu = mu_opt;
out.C = C_opt;
out.a = a_opt;
out.lambda = l_opt;
out.U = U_opt;
out.success = exitflag;

end