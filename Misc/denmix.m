function [x, psnr] = denmix(y, f, sigma, database)

%% Denoise using previously trained mixture

if ~iscell(database)
    load(database);
else

    if length(database) == 4
        prob = database{1};
        Scomp = database{2};
        Ucomp = database{3};
        supportvar = database{4};    
    elseif length(database) == 3
        prob = database{1};
        Scomp = database{2};
        Ucomp = database{3};
    end
end


pd = sqrt(size(Ucomp,1));
ps = 1;

% global_dc = mean(y(:));
global_dc = 0;
y = y - global_dc;
yy = wextend(2,'sym',y,[pd,pd]);
y_patches_ac = im2colstep(yy,[pd,pd],[ps,ps]);


if ~exist('supportvar', 'var')
%     fprintf('No fixed support.\n')
    supportvar = [];
end

[dimens,num] = size(y_patches_ac);
[m, n] = size(f);

y_patches_dc=mean(y_patches_ac);
y_patches_ac= bsxfun(@minus, y_patches_ac , y_patches_dc);

% 
max_im = max(max(y_patches_ac(:)), 0.5);
min_im = min(min(y_patches_ac(:)), -0.5);
scale = max_im - min_im;

y_patches_ac = (y_patches_ac)/scale;


[x_hat_patches] = GMM_inference(y_patches_ac,zeros(dimens,length(prob)),prob,Scomp,Ucomp,sigma/scale, 'supportvar', supportvar);

weights = ones(size(x_hat_patches));

x_hat_patches = min(max(x_hat_patches, -0.5 - global_dc),0.5-global_dc);
% 
x_hat_patches = (x_hat_patches)*scale;

x_hat_patches = bsxfun(@plus, x_hat_patches , y_patches_dc);

x_hat1 = col2imstep(x_hat_patches.*(weights),size(yy),[pd,pd],[ps,ps]);
normalize = col2imstep((weights),size(yy),[pd,pd],[ps,ps]);
x_hat1 = x_hat1 ./ normalize;

x_hat = x_hat1(pd+1:(m+pd),pd+1:(n+pd)) + global_dc;

x = x_hat;

x_hatb = x;
psnr = 10*log10(1/mean( (f(:)-x(:)).^2));

PSNR_best= psnr;

x = x_hatb;
psnr = PSNR_best;
end