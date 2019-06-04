function y = multinorm_svd(x,S,U)
% Evaluates a multidimensional Gaussian
% of mean m and covariance matrix C = U*S*U'
% at the array of points x
%
[dim, npoints] = size(x);
invS = 1./( S + eps);
in = U*diag(invS)*U';
%ff = ((2*pi)^(-dim/2))*((prod(S + realmin))^(-0.5));
a = prod(sqrt(invS));
ff = ((2*pi)^(-dim/2))* a;

var = 1e-5;
while all(isinf(a))
    invS = 1./( S + var);
    a = prod(sqrt(invS));
    ff = ((2*pi)^(-dim/2))* a;
    var = var*10;
end

%quadform = zeros(1,npoints);
%centered = (x-repmat(m,1,npoints));
if dim ~= 1
    y = ff * exp(-0.5*sum(x.*(in*x)));
%    y = ff * exp(-0.5*sum(bsxfun(@times,x,(in*x))));
    %y = ff * exp(-0.5*dot(x,(in*x)));
else
   y = ff * exp(-0.5*x.^2*in );
end

% 
% Sigma = U*diag(S)*U';
% [R,err] = cholcov(Sigma,0);
% 
% %// Create array of standardized data, and compute log(sqrt(det(Sigma)))
% xRinv = R\x;
% logSqrtDetSigma = sum(log(diag(R)));
% 
% %// Finally get the quadratic form and thus, the final output
% quadform = sum(xRinv.^2, 1);
% p_out = exp(-0.5*quadform - logSqrtDetSigma - dim*log(2*pi)/2);


