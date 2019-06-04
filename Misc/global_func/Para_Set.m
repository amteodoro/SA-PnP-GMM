function [name_image,band_remove,band_set,nr,nc,N_band,nb_sub,X_real,XH,XHd,XHd_int,XM,VXH,VXM,...
    psfY,psfZ_unk,s2y_real,s2z_real,SNR_HS,SNR_MS,miu_x_real,s2_real,P_inc,P_dec,eig_val]=Para_Set(seed,scale,subMeth,SNR_R, dataset, experiment)

%% Select the image source

if dataset == 1
    %nr=128;nc=128; name_image = 'moffet_ROI3_bis.mat';N_band=172;band_set=[26 17 10];nb=4; N_start_r=50;N_start_c=50;nb_sub=10;
    nr=128;nc=128;name_image = 'moffet_ROI3_bis.mat';N_band=176;band_set=[20 11 4];N_start_r=0; N_start_c=0;nb_sub=10;%For MS+HS ------
    % nr=64;nc=64;name_image = 'moffet_ROI3_bis.mat';N_band=3;N_start=0;nb_sub=3;band_set=[]; %For Pan+MS
    % nr=64;nc=64;name_image = 'moffet_ROI3_bis.mat';N_band=177;band_set=[20 11 4];N_start=0;nb_sub=4; %For Pan+HS
else
    nr=128;nc=128;name_image='pavia.mat';N_band=103-10;band_set=[45 25 8];N_start_r=50;N_start_c=50;nb_sub=5;  %For MS+HS -------
    % nr=512;nc=256;name_image='pavia.mat';N_band=103-10;band_set=[45 25 8];N_start_r=90;N_start_c=00;nb_sub=5;  %For MS+HS
    % nr=64;nc=64;name_image='pavia.mat';N_band=50;band_set=[45 25 8];N_start=50;nb_sub=5;%For Pan+HS
    % name_image='pleiades_subset.mat';N_band=103;
    % SNR_R=SNR_R_set(i_R);
    % SNR_R=inf;
end
%% Constructing the groundtruth image
[X_real,band_remove]= real_image(name_image,nr,nc,N_band,N_start_r,N_start_c);%X_temp=X_real;

%% Set the noise power of HS and MS data
if experiment == 1 || experiment == 3
    SNR_HS=[50*ones(N_band-50,1);50*ones(50,1)];
    
    SNR_MS=50;
else
    SNR_HS=[35*ones(N_band-50,1);30*ones(50,1)];
    
    SNR_MS=30;
end

if experiment == 1 || experiment == 2
    pan = 1;
else
    pan = 0;
end
%% Generate the HS and MS images
[psfY,psfZ_unk,XH,XM,s2y_real,s2z_real]= HS_MS(X_real,SNR_HS,SNR_MS,SNR_R,name_image,band_remove, pan);

%% Downsampled HS image
XHd=XH(1:psfY.ds_r:end,1:psfY.ds_r:end,:);
XHd_int=ima_interp_spline(XHd,psfY.ds_r);
%% HS subspace identification: Identifying the subspace where HS data live in
temp=reshape(XHd,[size(XHd,1)*size(XHd,2) N_band])';

if strcmp(subMeth,'Hysime')
    [w,Rn] = estNoise(temp, 'additive');
    [nb_sub,P_vec]=hysime(temp,w,Rn);
    eig_val = 1;
elseif strcmp(subMeth,'PCA')
    [P_vec,eig_val]=fac(XHd);
    PCA_ratio=sum(eig_val(1:nb_sub))/sum(eig_val); %[P_vec,eig_val]=fac(X_real);
    P_vec=P_vec(:,1:nb_sub); % Each column of P_vec is a eigenvector
end
if scale==1
    P_dec=diag(1./sqrt(eig_val(1:nb_sub)))*P_vec';
    P_inc=P_vec*diag(sqrt(eig_val(1:nb_sub)));
elseif scale==0
    P_dec=P_vec';
    P_inc=P_vec;
end
%% Project the image to subspace and project it back
miu_x_real=squeeze(mean(mean(X_real,1),2));
s2_real = cov(reshape(X_real,[size(X_real,1)*size(X_real,2) size(X_real,3)]));

VXH=reshape(XHd,size(XHd,1)*size(XHd,2),size(XHd,3))';
VXM=reshape(XM,size(XM,1)*size(XM,2),size(XM,3))';
