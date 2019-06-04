function [x_hat, global_post_cov, weight_factors] = GMM_inference(y,mu,prob,S,U,sigma,varargin)
%
% produces the MMSE estimate under the following
% observation model and prior:
%
% observation model:   y = x + n,
% where n is white Gaussian noise of variance sigma^2
%
% prior (Gaussian mixture):
%      p(x) = sum_{j=1}^K  prob(j) Normal(x|mu(j),cov(j))
%
%      where each covariance is represented by its SVD
%
% Inputs:
%   y (m by n): n points of dimension m, which are the
%               observations
%
%   prob (K by 1):  vector of mixing probabilities of the
%                   mixture prior
%
%   mu (m by K):  mean vector of each of the K components
%
%   S     (m by K): the vectors of singular values of each of
%                   the K covariance matrice.
%
%   U     (m by m by K): where U(:,:,j) contains the singular
%                        vectors of the j-th covariance matrix.
%
%   sigma:      the noise standard deviation
%
%   Optional input arguments
%
%   For scene-adaptation
%
%   supportvar = posterior weights from training
%
%   For block-GMMs
%
%   E = subspace vectors
%
%   pd = spatial dimension of the blocks, such that a block is of size
%           pd * pd * dimesion of subspace
%
% Outputs:
%   x_hat  (d by n): the estimates corresponding to each of thr
%                    observations in y.
%
%   global_post_cov (m by m): the posterior covariance.
%
%
%  Authors: Afonso Teodoro (afonso.teodoro91@gmail.com)
%           Mario A. T. Figueiredo  (mtf@lx.it.pt)


[d,K] = size(S);
[dimens,n] = size(y);

if numel(varargin)~=0
    for i=1:2:numel(varargin)
        if strcmp(varargin{i},'E')
            E = varargin{i+1};
        end
        if strcmp(varargin{i},'pd')
            pd = varargin{i+1};
        end
        if strcmp(varargin{i},'supportvar')
            supportvar = varargin{i+1};
        end
    end
end
if exist('pd') == 0
    pd = 1;
end
if exist('E') == 0
    E = eye(dimens);
end
if exist('supportvar') == 0
    supportvar = [];
end

d1 = size(E, 1);


sigma2 = sigma.^2;
if length(sigma) > 1
    eyesigma2 = kron(E*diag(sigma2)*E',eye(pd^2));
else
    eyesigma2 = kron((eye(d1)*sigma2)*E*E',eye(pd^2));
end
filters = zeros(d,d,K);

post_cov = zeros(d,K);
cov = zeros(size(U));

if ~isempty(supportvar)
    fixed_support = 1;
else
    fixed_support = 0;
end

for f=1:K

    cov(:,:,f) = U(:,:,f)*bsxfun(@times,S(:,f),U(:,:,f)');
    filters(:,:,f) = cov(:,:,f)/(eyesigma2 + cov(:,:,f));

    post_cov(:,f) = diag(cov(:,:,f) - filters(:,:,f)*cov(:,:,f));
end


if fixed_support == 1 % scene-adaptation: keep weights fixed
    weight_factors = supportvar;
else
    weight_factors = computePost(y, prob, mu, cov + repmat(eyesigma2,[1,1,K]));
    A = bsxfun(@minus,weight_factors,max(weight_factors,[],1));
    B = exp(A);
    weight_factors = bsxfun(@rdivide,B,sum(B,1));
end

x_hat = zeros(d,n);
x_hat_2 = zeros(d,n);

for j=1:K
    aux = bsxfun(@plus,(filters(:,:,j)*(y-repmat(mu(:,j),1,n))),mu(:,j));
    x_hat = x_hat + bsxfun(@times,aux,weight_factors(j,:));
    x_hat_2 = x_hat_2 + bsxfun(@times,bsxfun(@power,aux,2),weight_factors(j,:));
end

global_post_cov = post_cov * weight_factors + x_hat_2 - x_hat.^2;


