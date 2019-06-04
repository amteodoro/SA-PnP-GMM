function [bestpp,bestcovS,bestcovU,best_loglike_s, indic, supportvar] = ...
    EM_zeromean(y,k,sigma, varargin)
%
% y = data matrix  (number of rows = number of variables;
%                   number of colums = number of points)
%
% k = number of components
%
% sigma = noise standard deviation
%
% Optional input arguments for block-GMMs
%
% numim = number of bands in training
%
%
% bestpp = estimate of the component probabilities
% bestcovS = estimate of the covariances singular values
% bestcovU = estimate of the covariances singular vectors
% best_loglike_s = evolution of the loglikelihood in the best run
% supportvar = posterior weights; to be kept fixed in scene-adaptation
%
%
% Authors: Afonso Teodoro (afonso.teodoro91@gmail.com)
%          Mario A. T. Figueiredo  (mtf@lx.it.pt)
%
%


n=numel(varargin);
if n~=0
    for i=1:2:n
        if strcmp(varargin{i},'numim')
            numim = varargin{i+1};
        end
    end
end
if exist('numim') == 0
    numim = 1;
end

randn('seed', 0);

verb = 0;
[dimens,npoints] = size(y);
mult = npoints/numim;

ysub = y(:,1:numim:end);
sigma2 = sigma.^2;

regularization = eps*eye(dimens);

maxloglike = -realmax;
startn = 0;

indic = zeros(k,size(ysub,2));

cov = zeros(dimens,dimens,k);
randstarts = 1;
th = 1e-5;

while startn < randstarts
    
    iter = 1;
    loglike_s = [];
    startn = startn + 1;
    breakflag = 0;
    
    Ucomp = zeros(dimens,dimens,k);
    Scomp = zeros(dimens,k);
    prob = (1/k)*ones(k,1);
    globcov = (y*y')/npoints;
    %         rng('default')
    rand('seed', 0);
    
    for i=1:k
        estcovar = diag((0.5+0.3*rand(1,dimens))*0.95*max(diag(globcov)));
        [UUU,SSS] = svd(estcovar);
        Ucomp(:,:,i) = UUU;
        Scomp(:,i) = diag(SSS);
    end
    
    
    
    maxloglike = -realmax;
    
    prevloglike = -realmax;
    cont = 1;
    while(cont) % EM loop
        
        for mod=1:k
            indic(mod,:) = prob(mod)*multinorm_svd_zeromean(ysub,...
                Scomp(:,mod)+sigma2,Ucomp(:,:,mod));
        end
        
        normindic = bsxfun(@times,indic,1./(realmin + sum(indic,1)));
        
        
        for mod=1:k
            
            
            aux = bsxfun(@times,ysub,normindic(mod,:)/(eps+sum(normindic(mod,:))));
            estcov = aux*ysub';
            
            cov(:,:,mod) = estcov;
            [Ucomp(:,:,mod),S, ~] = svd(cov(:,:,mod)+regularization);
            Scomp(:,mod) = max(diag(S)-sigma2,0);
            
        end
        
        prob = sum(normindic,2) / npoints;
        prob = prob/sum(prob);
        
        
        loglike = sum(log(sum(indic)));
        
        deltlike = loglike - prevloglike;
        if (verb~=0) %&& mod(iter,5) == 0
            fprintf(...
                'iter = %g, k = %g, logL = %0.5g, deltalogL/th = %0.7g\n', ...
                iter, k, loglike, abs(deltlike/loglike)/th);
        end
        
        if abs(deltlike/loglike) < th
            cont = 0;
        end
        
        loglike_s(iter) = loglike;
        iter = iter + 1;
        
        if iter > 30
            cont = 0;
        end
        
        if iter > 10 && isnan(abs(deltlike/loglike)/th)
            cont = 0;
            startn = startn - 1;
            loglike = -realmax;
            prob = 0;
            sigma = sigma*2;
            sigma2 = sigma.^2;
        end
        
        if breakflag == 1
            cont = 0;
            loglike = -realmax;
            startn = startn - 1;
        end
        
        prevloglike = loglike;
        
    end % the EM loop
    
    if loglike > maxloglike
        best_loglike_s = loglike_s;
        maxloglike = loglike;
        bestcovS = Scomp;
        bestcovU = Ucomp;
        bestpp = prob;
    end
    
    supportvar2 = zeros(length(bestpp),npoints);
    
    for f=1:length(bestpp)
        
        supportvar2(f,:) = bestpp(f)*multinorm_svd_zeromean(y,bestcovS(:,f)+sigma2,bestcovU(:,:,f));
    end
    supportvar2 = bsxfun(@times,supportvar2,1./(realmin + sum(supportvar2,1)));
    supportvar = reshape(supportvar2, [length(bestpp), mult, numim]);
    supportvar = sum(supportvar, 3)/numim;
    
    
end % the random starts loop



