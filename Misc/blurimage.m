function [f_blur, R, RT, sigma, hshift, h, f_clean] = blurimage(input, experiment_number, varargin)

[fM,fN]=size(input);


if nargin == 3
    experiment_number = 1000;
end

switch experiment_number
    case 1
        sigma=sqrt(2)/255;
        for x1=-7:7; for x2=-7:7; h(x1+8,x2+8)=1/(x1^2+x2^2+1); end, end; h=h./sum(h(:));
    case 2
        sigma=sqrt(8)/255;
        s1=0; for a1=-7:7; s1=s1+1; s2=0; for a2=-7:7; s2=s2+1; h(s1,s2)=1/(a1^2+a2^2+1); end, end;  h=h./sum(h(:));
    case 3
        BSNR=40;
        sigma=-1; % if "sigma=-1", then the value of sigma depends on the BSNR
        for x1=-4:4; for x2=-4:4; h(x1+5,x2+5)=1; end, end; h=h./sum(h(:));
    case 4
        sigma=7/255;
        h(1:5,1:5)=[1 4 6 4 1]'*[1 4 6 4 1]; h=h./sum(h(:));  % PSF
    case 5
        sigma=2/255;
        h(1:25,1:25)=fspecial('gaussian', 25, 1.6);
        %h = fftshift(h);
    case 6
        sigma=8/255;
        h(1:25,1:25)=fspecial('gaussian', 25, .4);
        %h = fftshift(h);
        %extra experiments
    case 7
        BSNR=30;
        sigma=-1;
        h=ones(9); h=h./sum(h(:));
    case 8
        BSNR=20;
        sigma=-1;
        h=ones(9); h=h./sum(h(:));
    case 9
        BSNR=40;
        sigma=-1;
        h=fspecial('gaussian', 25, 1.6);
    case 10
        BSNR=20;
        sigma=-1;
        h=fspecial('gaussian', 25, 1.6);
    case 11
        BSNR=15;
        sigma=-1;
        h=fspecial('gaussian', 25, 1.6);
    case 12
        BSNR=40;
        sigma=-1; % if "sigma=-1", then the value of sigma depends on the BSNR
        h=ones(19); h=h./sum(h(:));
    case 13
        BSNR=25;
        sigma=-1; % if "sigma=-1", then the value of sigma depends on the BSNR
        h=ones(19); h=h./sum(h(:));
    otherwise
        sigma = 0;
        %         BSNR = 40;
        h = varargin{1};
end

hpad = zeros(fM,fN);
hpad(1:size(h,1),1:size(h,2)) = h;

hshift = circshift(hpad, [-floor(size(h,1)/2), -floor(size(h,2)/2)]);


% define the function handles that compute
% the blur and the conjugate blur.
R = @(x) real(ifft2(fft2(hshift).*fft2(x)));
RT = @(x) real(ifft2(conj(fft2(hshift)).*fft2(x)));

f_blur = R(input);


if sigma == -1;   %% check whether to use BSNR in order to define value of sigma
    sigma=sqrt(norm(f_blur(:)-mean(f_blur(:)),2)^2 /(fM*fN*10^(BSNR/10)));
    %     Xv% compute sigma from the desired BSNR
end
