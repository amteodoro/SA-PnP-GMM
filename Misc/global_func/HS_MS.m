%Simulate the Hyperspectral image and Multispectral image
function [psfY,psfZ,Y,Z,sigma2y_real,sigma2z_real] = HS_MS(X_real,SNR_HS,SNR_MS,SNR_R,name_image,band_remove, pan)
randn('seed',0)

L=size(X_real,3);
if strcmp(name_image,'moffet_ROI3_bis.mat')
%     psfY = fspecial('average',[4 4]);
    % psf = fspecial('gaussian',[5 5],1.5);
    % psf = 1+abs(fspecial('sobel')); 
    if L<=10
        psfZ_real = ones(1,L)/L;
        %     psfZ = [1/6 1/2 1/3];
    elseif L>10
        %     psfZ = kron(eye(10),ones(1,L/10)/L*10);
%         psfZ_real = ones(1,L);
% %         psfZ_real = [ones(1,3)/3 zeros(1,L-3)];
%         psfZ_real=psfZ_real./repmat(sum(psfZ_real,2),1,size(psfZ_real,2));
        
        psfZ_real = ones(1,L)/L;
       load('panF2full.mat')
       [no_wa, ~] = size(out);
       xx  = linspace(1,no_wa,180);
       x = 1:no_wa;
       outcrop = spline(x, out(:,2), xx);
            load('F2.mat')

%         Construct the spectral response for LANDSAT image
        R=zeros(5,224);
        R(2,7:14)=mean(reshape(rep6S(16:47,3),[4 8])); %% Averaging every 4 bands
        R(3,13:24)=mean(reshape(rep6S(41:88,4),[4 12]));
        R(4,24:32)=mean(reshape(rep6S(86:121,5),[4 9]));
        R(5,40:57)=mean(reshape(rep6S(136:207,6),[4 18]));
        R(1,13:57)=mean(reshape(outcrop(:),[4 45])); % PAN
%         psfZ_real(5,122:149)=mean(reshape(rep6S(445:556,7),[4 28]));
%         psfZ_real(6,175:213)=mean(reshape(rep6S(641:796,8),[4 39]));
        R(:,band_remove)=[];
%         psfZ_real(7,2:51)=0.98+rand(1,50)*0.01;
        if pan == 1
            psfZ_real = R(1,:);
        else
            psfZ_real = R(2:5,:);
        end
        psfZ_real=psfZ_real./repmat(sum(psfZ_real,2),1,size(psfZ_real,2));
    end
elseif strcmp(name_image,'pavia.mat')
%     load('pavia.mat');clear B Yh Ym Z mask;
%     Btemp = fftshift(B);
%     psfY.B=Btemp(302:310,165:173);
    % USE IKNOS RESPONSES
    % spectral responses (wavelenghths, pan, blue, green, red, NIR, in nanometers)
    load iknos_spec_resp.mat
    [no_wa, ~] = size(iknos_sp);

    % map iknos wavelengths  into rosis (430 - 860 nm) bands
    % find valid interval iknos \subset rosis
    valid_ik_bands = [];
    for i = 1:no_wa
        if iknos_sp(i,1) >=  430 &&  iknos_sp(i,1) <= 860
            valid_ik_bands = [valid_ik_bands i];
        end
    end
    no_wa = length(valid_ik_bands);

    % spline interpolation
    xx  = linspace(1,no_wa,L+10);
    x = 1:no_wa;

    % 1 pan; 2 - blue; 3 - green; 4 - red; 5 - NIR
    for i=1:5
        R(i,:) = spline(x,iknos_sp(valid_ik_bands,i+1),xx);
    end
    
%     load('R.mat');        
        if pan == 1
            psfZ_real=R(1,1:end-10);
        else
            psfZ_real=R(2:5,1:end-10);
        end
    
%     psfZ_real=psfZ_real(1:3,:);
    psfZ_real=psfZ_real./repmat(sum(psfZ_real,2),1,size(psfZ_real,2));

%      psfZ_real = ones(1,L)/L;

%     temp=size(R,2)-mod(size(R,2),4);
%     psfZ_real=zeros(4,size(R,2));
%     for i=1:4
%         psfZ_real(i,(i-1)*temp/4+1:i*temp/4)=ones(1,temp/4);
%     end
%     psfZ_real(4,temp+1:end)=1;
%     psfZ_real=psfZ_real./repmat(sum(psfZ_real,2),1,size(psfZ_real,2));

%    
elseif strcmp(name_image,'pleiades_subset.mat')
    load pleiades_subset.mat;
%   psfZ = ones(1,L)/L;
    psfZ_real=[0.077919      0.37514      0.50553    0.0057739]; 
end
%% Generate the downsampling matrix 'mask' and blurring mask 'B'
% load B_mask
% B = circshift(B, [20 20]);
% B = B(1:size(X_real,1),1:size(X_real,2));
% psfY.B= circshift(B, [-20 -20]);

% blurr matrix Type 1: average filter; Type 2: expotienal degrading
blur_type=2;
psfY.B=gen_degr_mat(blur_type,size(X_real,1),size(X_real,2),2,2);
% psfY.dsp=mask(1:size(X_real,1),1:size(X_real,2));
psfY.ds_r=4;
mask=zeros(size(psfY.B));
mask(1:psfY.ds_r:end,1:psfY.ds_r:end,:)=1;
psfY.dsp=mask;
%% Spatial blurring
Y = func_blurringY(X_real,psfY);       
% Y_tem=Y(1:psfY.ds_r:end,1:psfY.ds_r:end,:);Ps = mean(Y_tem(:).^2);
Ps = squeeze(mean(mean(Y(1:psfY.ds_r:end,1:psfY.ds_r:end,:).^2,1),2));
sigma2y_real = Ps.*(10.^(-SNR_HS/10));  %Caclulate the noise power
% sigma2y_real = (sum(Y(:).^2)/(10^(SNR_HS/10)*numel(Y)));
% rng(100,'v5normal'); %Set the seed
for i=1:L
    Y(:,:,i) = Y(:,:,i) + randn(size(Y(:,:,i)))*sqrt(sigma2y_real(i));
end
% Y = Y + randn(size(Y))*sqrt(sigma2y_real);
% dsp_op=zeros(size(Y));
% dsp_op(1:psfY.ds_r:end,1:psfY.ds_r:end,:)=1;
Y = Y.*repmat(psfY.dsp,[1 1 size(Y,3)]);

% %% Unregisration
% n_off=8;
% Y=[Y(:,n_off+1:end,:) Y(:,1:n_off,:)];

%% Spectral blurring: psfZ is the real observation matrix while psfZ_unk is contamined with noise
% rng(130,'v5normal'); %Set the seed
% temp=psfZ; temp(psfZ>0)=psfZ_tem(psfZ>0); psfZ_unk=temp./kron(sum(temp),ones(size(temp,1),1));
sig_F=sqrt(10^(-SNR_R/10)*norm(psfZ_real,'fro')^2/numel(psfZ_real));
psfZ=psfZ_real+sig_F*randn(size(psfZ_real)).*(psfZ_real~=0);

Z = func_blurringZ(X_real,psfZ_real);
% Ps = mean(Z(:).^2)*ones(size(Z,3),1);
Ps = squeeze(mean(mean(Z.^2,1),2));
sigma2z_real = Ps.*(10.^(-SNR_MS/10));
% sigma2z_real = (sum(Z(:).^2)/(10^(SNR_MS/10)*numel(Z)));
% rng(120,'v5normal'); %Set the seed
for i=1:size(Z,3)
    Z(:,:,i) = Z(:,:,i) + randn(size(Z(:,:,i)))*sqrt(sigma2z_real(i));
end
% Z = Z + randn(size(Z))*sqrt(sigma2z_real);

for i=1:size(Z,3)
    Im = Z(:,:,i);
    Zaff(:,:,i) = Z(:,:,i)/sqrt(mean(Im(:).^2));
end
% Zaff = abs(Z/max(Z(:)));
% figure(1)
% subplot(1,3,1);imshow(X_real);   title('Reference')
% subplot(1,3,2);imshow(Yaff);title('Hyperspectral Image')
% subplot(1,3,3);imshow(Zaff);title('Multispectral Image')