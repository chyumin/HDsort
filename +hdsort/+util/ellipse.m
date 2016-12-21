function M = ellipse(mu, Cov_, varargin)

    [U lambda] = eig(Cov_);
    %fit.lambda_pos(1), fit.lambda_pos(2), sign(diff(fit.lambda_pos))*atan(fit.U(2,1)/fit.U(1,1)));
    lambda = sort(diag(lambda), 'descend');
    assert( ~any( lambda <= 0 ), 'Eigenvalues must be strictly positive!')
    a = sqrt(lambda(1)); b = sqrt(lambda(2));
    alpha_ = atan(U(1,1)/U(2,1));
 
    th = 0:pi/50:2*pi;    
    M_ = [cos(th)*a; sin(th)*b]';
    R = [cos(alpha_) -sin(alpha_); sin(alpha_) cos(alpha_)];
    M = M_*R + repmat([mu(1) mu(2)], [size(M_,1) 1]);
    
end