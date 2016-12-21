function [out P] = receptiveFieldSVD(RF, varargin)
% Input:
%  - RF: receptive field as 2D-matrix.   
%  - varargin: [so far only 'debug']
%
% Output:
%  - out: structure containing all fitting parameters
%  - P: input preferences

P.debug = false;
%P.MaxIter = [];
P = mysort.hdsort.util.parseInputs(P, varargin);

xrange = 1:size(RF,1);
yrange = 1:size(RF,2);

% out.c = median(RF(:));
% RF = RF - out.c;
% 
% if P.debug
%     figure; imagesc(RF);
%     figure;
% end
% 
% EH = @(params) E(params, RF, xrange, yrange);
% 
% [max_x max_y] = mysort.hdsort.util.matrixArgMax(abs(RF));
% ampmax = RF(max_x,max_y);
% params0 = [max_x max_y 10*rand(1) 10*rand(1) 10*rand(1) 0.5 ampmax mean(RF(:))];
% 
% if ~isempty(P.MaxIter)
%     [paramsopt, fopt, exitflag] = fminsearch(EH, params0, optimset('MaxIter', P.MaxIter));
% else
%     [paramsopt, fopt, exitflag] = fminsearch(EH, params0); 
% end
% 
% 
%     function [mu Cpos Cneg a b U l_pos l_neg ] = interpretParameters(params)
%         mu = params(1:2);
%         %l_pos = abs(params(3:4));
%         
%         %% Constrain the difference between l_pos(1) and l_pos(2) to a ration of epsilon:
%         if 1
%             epsilon = 0.1;
%             lower_bound = log(epsilon)/log(10);
%             upper_bound = abs(lower_bound);
%             
%             l_pos(1) = abs(params(3));
%             l_pos(2) =  10^( (upper_bound-lower_bound)/( 1+exp(- params(4) ) ) + lower_bound) * l_pos(1) ;
%         else
%             epsilon = 0.1;
%             if l_pos(1)/l_pos(2) < epsilon
%                 l_pos(1) = l_pos(2)*epsilon;
%             elseif l_pos(2)/l_pos(1) < epsilon
%                 l_pos(2) = l_pos(1)*epsilon;
%             end
%         end
%         
%         l_neg = abs(params(5)) * l_pos;
%         
%         u11 = sin(params(6));
%         u22 = sqrt(1 - u11^2);
%         u1 = [u11 u22];
%         u2 = [-u1(2) u1(1)];
%         
%         U = [u1; u2];
%         Lpos = [l_pos(1) 0; 0 l_pos(2)];
%         Cpos = U*Lpos*U';
%         
%         
%         Lneg = [l_neg(1) 0; 0 l_neg(2)];
%         Cneg = U*Lneg*U';
%         
%         assert( ~any([Cpos(1,1), Cpos(2,2), Cneg(1,1), Cneg(2,2)] < 0), 'Cov matrices not well defined!')
%          
%         a = params(7);
%         b = -sign(a)*abs(params(8));
%     end
% 
%     function e = E(params, RF, xrange, yrange)
%         
%         [mu Cpos Cneg a b ] = interpretParameters(params);
%         
%         RF_fit = RFmodel(xrange, yrange, mu, Cpos, Cneg, a, b);
%         
%         e = norm(RF_fit(:) - RF(:));
%         
%         if P.debug
%             e
%             subhdsort.plot.2,1,1)
%             %imagesc(reshape(RF_fit(:), size(RF,1), size(RF,2)) )
%             imagesc(RF_fit)
%             subhdsort.plot.2,1,2)
%             %imagesc(reshape(RF(:), size(RF,1), size(RF,2)) )
%             imagesc(RF)
%             pause(0.1)
%         end
%     end
% 
%     function rf = RFmodel(xrange, yrange, mu, C1, C2, a, b)
%         [X Y] = meshgrid(xrange, yrange);
%         XX = [X(:) Y(:)];
%         rf = mvnpdf(XX, mu, C1)*a - mvnpdf(XX, mu, C2)*b;
%         rf = reshape(rf, length(yrange), length(xrange))';
%     end
% 
% [mu_opt Cpos_opt Cneg_opt a_opt b_opt U_opt l_pos_opt l_neg_opt ] = interpretParameters(paramsopt);
% %out.RF = reshape(RFmodel(xrange, yrange, mu_opt, Cpos_opt, Cneg_opt, a_opt, b_opt), size(RF,2), size(RF,1))';
% out.RF = RFmodel(xrange, yrange, mu_opt, Cpos_opt, Cneg_opt, a_opt, b_opt);
% 
% if P.debug
%     figure;
%     subhdsort.plot.2,1,1); imagesc(RF);
%     subhdsort.plot.2,1,2); imagesc(out.RF)
% end
% 
% out.mu = mu_opt;
% out.Cpos = Cpos_opt; out.Cneg = Cneg_opt;
% out.a = a_opt; out.b = b_opt;
% out.lambda_pos = l_pos_opt;
% out.lambda_neg = l_neg_opt;
% out.U = U_opt;
% out.success = exitflag;

end