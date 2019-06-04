%% Name: demo_SA_PnP_GMM
%
%  Code used to reproduce the fusion experiments presented in
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

close all
clear

dataset = 1;        % 1 - moffet; 2 - pavia

experiment = 1;     % 1 - Train GMM from PAN, 50dB
                    % 2 - Train GMM from PAN, 30dB
                    % 3 - Train GMM from 4 MS bands, 50dB
                    % 4 - Train GMM from 4 MS bands, 30dB

% define random states
rand('state',10);
randn('state',10);

% ADMM iterations
iters = 100;

generate=1;subMeth='PCA';FusMeth='Sparse';
scale=1;
seed=1;

switch dataset
    
    case 1
        if experiment == 1
            database = 'moffet_ROI3_50db_pan';
        elseif experiment == 2
            database = 'moffet_ROI3_35db_pan';
        elseif experiment == 3
            database = 'moffet_ROI3_50db_ms';
        else
            database = 'moffet_ROI3_35db_ms';
        end
        loaddataset = 1;
    case 2
        if experiment == 1
            database = 'pavia_50db_pan';
        elseif experiment == 2
            database = 'pavia_35db_pan';
        elseif experiment == 3
            database = 'pavia_50db_ms';
        else
            database = 'pavia_35db_ms';
        end
        loaddataset = 2;
end

d = 4; % decimation factor
meanval = 0;
maxval = 1;
minval = 0;
PSNR=[];

%% Generate the data
[name_image,band_remove,band_set,nl,nc,L,p1,Zim,Yhim,Yhd,Yhd_int,Ymim,VXH,VXM,psfY,psfZ_unk,...
    sigmah2,sigmam2,SNR_HS,SNR_MS,miu_x_real,s2_real,P_inc,P_dec,eig_val]=Para_Set(seed,1,subMeth,inf,loaddataset,experiment);

PAN = Ymim;
np = nl*nc;
Yh = reshape(Yhim,size(Yhim,1)*size(Yhim,2),L)';

E = P_inc;

p = size(E,2)

pinvE = (E'*E)\E';

Z = reshape(Zim,nl*nc,L)';
sigmah2m = mean(sigmah2);
sigmam2m = mean(sigmam2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define a circular convolution (the same for all bands) accepting a
% matrix  and returnig a matrix

ConvC = @(X,FK)  reshape(real(ifft2(fft2(reshape(X', nl,nc,p)).*repmat(FK,[1,1,p]))), nl*nc,p)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert matrix to image
conv2im  = @(X)  reshape(X',nl,nc,p);
% convert image to matrix
conv2mat = @(X)  reshape(X,nl*nc,p)';

            
% define horiontal difference operator kernel
dh = zeros(nl,nc);
dh(1,1) = 1;
dh(1,nc) = -1;

dv = zeros(nl,nc);
dv(1,1) = 1;
dv(nl,1) = -1;

FDH = fft2(dh);
FDHC = conj(FDH);
FDV = fft2(dv);
FDVC = conj(FDV);

FB=fft2(psfY.B);

mask = zeros(size(psfY.B));
mask(1:d:end,1:d:end,:)=1;
dsp=mask;
mask = repmat(dsp(:)', p,1);

FBC = conj(FB);

IBD_B  = FBC  ./(abs(FB.^2) + 2);
IBD_II = 1    ./(abs(FB.^2) + 2);

R = psfZ_unk;

Ym = reshape(Ymim,nl*nc,size(Ymim,3))';


% ADMM parameters
lambda_phi = 1e-5; %[0.000001 0.00001 0.0001];

lambda_m = 1; %logspace(-2, 1, 4);

mu = 1e-3; %logspace(-3, 0, 4);


% Train a mixture if there isn't one saved
if exist(strcat(database, '.mat'), 'file')
    load(database);
    fprintf('Loading mixture...\n')
    if exist('supportvar', 'var')
        database = {prob, Scomp, Ucomp, supportvar};
    else
        database = {prob, Scomp, Ucomp};
    end
else
    fprintf('Training mixture...\n')
    pd = 8;
    K = 20;
    ps = 1;
    y_patches= [];
    clear xx
    for j = 1:size(PAN,3)
        xx(:,:,j) = wextend(2,'sym',PAN(:,:,j),[pd,pd]);
        y_patches = [y_patches, im2colstep(xx(:,:,j),[pd,pd],[ps,ps])];
    end
    
    y_patches_dc=mean(y_patches);
    y_patches_ac= bsxfun(@minus, y_patches , y_patches_dc);
    
    scale = 1;
    [prob,Scomp,Ucomp,~, ~, supportvar] = ...
        EM_zeromean(y_patches_ac,K,...
        sqrt(sigmam2m)/scale, 'numim', size(PAN,3));
    
    save(database, 'prob', 'Scomp', 'Ucomp', 'supportvar')
end

ETE = E'*E;
invETE = inv(ETE);

randn('seed',0)

fprintf('Parameters: mu = %f;\t lambda_m = %f;\t lambda_phi = %f;\n', mu, lambda_m, lambda_phi)

% auxiliary matrices
IE = inv(E'*E+mu*eye(p));
yyh = E'*Yh;

IRE = inv(E'*R'*R*E+mu/lambda_m*eye(p));
yym = E'*R'*Ym;

% define and initialize variables

U = pinvE*reshape(Yhd_int,nl*nc,size(Yhd_int,3))';%zeros(nl*nc,p)';%
V1 = ConvC(U,FB);
D1 = zeros(nl*nc,p)';
V2 = U;
D2 = D1;
V3 = U;
D3 = D1;
  
V3im = conv2im(V3);

reverseStr = '';

for i=1:iters
    %   min   ||UB - V1 - D1||_F^2  +
    %    U    ||U  - V2 - D2||_F^2  +
    %         ||U  - V3 - D3||_F^2
    
    U =    ConvC(V1+D1, IBD_B) + ConvC(V2+D2, IBD_II) + ConvC(V3+D3, IBD_II);
    
    %  max (1/2)||Y-EV1M|_F^2 + (mu/2)||UB - V1 - D1||_F^2
    %   V1
    NU1 =  ConvC(U,FB) - D1;
    V1 = IE*(yyh + mu*NU1).*mask + NU1.*(1-mask);
    
    %  max (lambda_m/2)||Y-REV2|_F^2 + (mu/2)||U - V2 - D2||_F^2
    %   V1
    NU2 =  U - D2;
    V2 = IRE*(yym + mu/lambda_m*NU2);
    
    % min lambda_phi phi(V3)+ (mu/2)||U - V3 - D3||_F^2
    % V3
    NU3 =  U-D3;
    V3 = NU3;
    
    V3im = conv2im(NU3);
    for k=1:p % Plugged-in GMM-denoiser on each coefficient band
        auxim = V3im(:,:,k);
        
        noiseVal(k) = sqrt(lambda_phi/mu);
        
        [auxim] = denmix(auxim, auxim, noiseVal(k), database);
        
        V3im(:,:,k) = auxim;
    end
    V3 = conv2mat(V3im);
    
    if mod(i, 10) == 0 || i == 1
        
        Uim = conv2im(U);
        num = sum(Zim .* reshape((E*U)',nl,nc,L), 3);
        den = sqrt(sum(Zim.^2, 3) .* sum(reshape((E*U)',nl,nc,L).^2, 3));
        SAM = sum(sum(acosd(num ./ den)))/(size(Zim,1)*size(Zim,2));
        
        Zhat = E*U;
        Zimhat = reshape(Zhat',nl,nc,L);
        Zim1 = Zim*(maxval+minval)-meanval;
        Zimhat1 = Zimhat*(maxval+minval)-meanval;
        aux = Zimhat1 - Zim1;
        
        for j = 1:L
            auxi = aux(:,:,j);
            auxZi = Zim1(:,:,j);
            PSNR(j) = 10*log10(max(auxZi(:))^2/(sum(auxi(:).^2)/(nl*nc)));
        end
        psnrrun = mean(PSNR);
        
        msg = sprintf('iter = %d, ||UB-V1|| = %2.2f, ||U-V2|| = %2.2f, ||UDH-V3|| = %2.2f; SAM = %4.2f; PSNR = %4.2f\n', ...
            i,  norm(NU1+D1-V1, 'fro'), norm(NU2+D2-V2, 'fro'), ...
            norm(NU3+D3-V3,'fro'), SAM, psnrrun);
        fprintf([reverseStr, msg]);
        reverseStr = '';
    end
    
    if mod(i, 1) == 0 || i == 1
        msg = sprintf('Iter = %d\n', i);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    % update Lagrange multipliers
    D1 = -NU1 + V1;    % D1 -(UB-V1)
    D2 = -NU2 + V2;    % D2 -(U-V2)
    D3 = -NU3 + V3;    % D3 -(U-V3)
    
end

Zhat = E*U;
Zimhat = reshape(Zhat',nl,nc,L);

% Denormalization
Zim1 = Zim*(maxval+minval)-meanval;
Zimhat = Zimhat*(maxval+minval)-meanval;
Z1 = reshape(Zim1,nl*nc,L)';

aux = Zimhat - Zim1;

Out = QualityIndices(Zimhat,Zim1,d);

for i = 1:L
    auxi = aux(:,:,i);
    auxZi = Zim1(:,:,i);
    PSNR(i) = 10*log10(max(auxZi(:))^2/(sum(auxi(:).^2)/(nl*nc)));
end

psnrrun = mean(PSNR)

resultsAll = [mu, lambda_m, lambda_phi, psnrrun, Out.sam, Out.ergas];


