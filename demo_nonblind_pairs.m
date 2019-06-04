%% Name: demo_nonbling_pairs
%
%  Code used to reproduce the noisy/blurred fusion experiments presented in
%
%  A. M. Teodoro, J. M. Bioucas-Dias and M. A. T. Figueiredo, 
%   "A Convergent Image Fusion Algorithm Using Scene-Adapted 
%   Gaussian-Mixture-Based Denoising," in IEEE Transactions 
%   on Image Processing, vol. 28, no. 1, pp. 451-463, Jan. 2019.
%
%
%  Author: Afonso Teodoro (afonso.teodoro91@gmail.com)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
close all

% Blur kernel
experiment_number = 3;

x=im2double(imread('Cameraman256.png'));
orig = x;


[blurred, ~, ~, sigmab, hshift, h] = blurimage(orig, experiment_number);

randn('seed',0)

blurred = blurred + sigmab*randn(size(blurred));

filter = h;

psnr_y = PSNR(orig,blurred,1,0);

sigma = [25];

noisy = orig + (sigma/255)*randn(size(orig));
sigma_hat = sigma/255;


pd = 8;
K = 20;
ps = 1;

% Training step using noisy image

yy = wextend(2,'sym',noisy,[pd,pd]);

patches = im2colstep(yy,[pd,pd],[ps,ps]);

patches_dc=mean(patches);
patches= bsxfun(@minus, patches , patches_dc);

[prob,Scomp,Ucomp, ~,~, supportvar] = ...
    EM_zeromean(patches,K,sigma_hat);

database = {prob; Scomp; Ucomp; supportvar};
filename = 'mix';
save(filename, 'prob', 'Scomp', 'Ucomp', 'supportvar');

database = filename;
isnr = [];

[ze, PSNR3] = denoising('clean', orig, 'noisy', noisy, 'sigma', sigma, 'sub', 0, 'dis', 0, 'external', database);

% Assume known blurring kernel
[~, R, RT, ~, hshift1, h1] = blurimage(orig, 90, h);

rho = 0.01;%[0.001 0.01 0.1 1 10];
lambda = 0.0001;%[0.00001 0.0001 0.001 0.01 0.1];
tau = 1;%[0.000001 0.00001 0.0001 0.001];

v3 = ze;

psnrv3 = [];
psnrx = [];
psnrv2 = [];

best_psnr = -inf;
fprintf('Parameters - rho: %f; \t lambda: %f; \t tau: %f;\n', rho, lambda, tau)

H_FFT = fft2(hshift1);
H2 = abs(H_FFT).^2;
filter_FFT = H2./(H2 + (rho) + lambda);
invLS = @(x, mu, filter_FFT) (1/mu)*( x - real( ifft2( filter_FFT.*fft2( x ) ) ) );

iter = 50;
x0 =  zeros(size(ze));
for i = 1:iter
    if i == 1 % Initializations
        x = x0;
        v1 = x;
        v2 = x;
        v3 = x;
        d1 = 0;
        d2 = 0;
        d3 = 0;
        y1 = noisy;
        y2 = blurred;
        RTy = RT(y2);
    end
    
    r = RTy + rho*(v3 + d3) + lambda*y1;
    x = invLS(r, rho + lambda, filter_FFT);
    
    auxim = x-d3;
    sigma_hat = NoiseEstimation(auxim, 8);
    if sigma_hat < 2/255
        sigma_hat = 2/255;
    end
    
    [v3, ~] = denmix(auxim, orig, sigma_hat, database);
    
    d3 = d3 - (x - v3);
    psnrv3(i) = PSNR(orig, v3, 1, 0);
    
    if psnrv3(i) > best_psnr
        best_psnr = psnrv3(i);
        best_i = i;
    end
    
    if ~mod(i,5) || i == 1
        psnrx(i) = PSNR(orig, x, 1, 0);
        psnrv2(i) = 1;%PSNR(orig, v2, 1, 0);
        fprintf('Iteration number: %d; \t PSNR v3: %4.2f; \t PSNR x: %4.2f; \t PSNR v2: %4.2f; \t ISNR: %4.2f\n', i, psnrv3( i), psnrx(i), psnrv2(i), psnrv3(i) - psnr_y)
    end
end



