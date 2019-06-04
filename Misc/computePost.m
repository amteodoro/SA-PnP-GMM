function [R, T] = computePost(X, w, mu, Sigma)

n = size(X,2);
k = size(mu,2);
R = zeros(k,n);
for i = 1:k
    R(i,:) = loggausspdf(X,mu(:,i),Sigma(:,:,i));
end
R = bsxfun(@plus,R,log(w));
T = logsumexp(R,2);
% R = exp(bsxfun(@minus,R,T));