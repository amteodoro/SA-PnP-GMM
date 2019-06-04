function [ze, PSNRbest, mssim] = denoising(varargin)

%% Default Parameters
seednum = 0;
ro = 0.03;
pdref = 8;
Kref = 20;
prob = 0;
Scomp = 0;
Ucomp = 0;
PSNRbest = zeros(1,2);
mode = 1; % use true noise variance
class = 0;
sub = 0;
% Step between patches
ps = 1;

cont = 1;
dis = 0; % Sharing the disagreement

mixture = [];
noisy = [];

if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'SUB'
                sub = varargin{i+1};
            case 'NOISY'
                noisy = varargin{i+1};
            case 'CLEAN'
                x = varargin{i+1};
            case 'SEED'
                seednum = varargin{i+1};
            case 'PD'
                pdref = varargin{i+1};
            case 'PS'
                ps = varargin{i+1};
            case 'K'
                Kref = varargin{i+1};
            case 'DIS'
                dis = varargin{i+1};
            case 'CONT'
                cont = varargin{i+1};
            case 'EXTERNAL'
                mixture = varargin{i+1};
                load(mixture)
                l = 1;
                for j = 1:4
                    p{j} = prob;
                    S{j} = Scomp;
                    U{j} = Ucomp;
                end
                pdref = sqrt(size(Ucomp,1));
                Kref = length(prob);
            case 'CLASS'
                class = 1;
            case 'SIGMA'
                sigma = varargin{i+1};
                %                 if sigma < 1
                %                     sigma = sigma*255;
                %                 end
                mode = 0;
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end
if ~exist('supportvar', 'var')
    %     fprintf('No fixed support.\n')
    supportvar = [];
end
randn('seed', seednum);

%%
iter = 1;
while cont == 1 % Oracle iterations, with signal boosting
    cont = 1;
    if isempty(noisy)
        y = x + sigma/255*randn(size(x));
    else
        y = noisy;
    end
    
    if iter ~= 1
        y = (1-ro)*y + ro*prevze;
    end
    
    %  PSNR computations
    PSNR_noisy    = -10*log10(var( (x(:)-y(:))));
        
    if sub
        [ze_up] = subsampDenoising(y, x, sigma, Kref); % Denoise lower resolution version of the image
        
        [sigma_hat] = noise_estimation(y, mode, sigma, 4, 1); % Estimate noise level; mode = 1 uses the true noise variance
        
        [ze, psnrOut] = denfun_global(x, y, sigma_hat, pdref, ps, Kref, mixture, dis, class, supportvar); % Denoise input image
                
        for i = 1:size(ze,3)
            ysub(:,:,i) = subsamp(ze(:,:,i));
            
            yzdif(:,:,i) = ze(:,:,i) - imresize(upsamp(ysub(:,:,i)), size(y(:,:,i)), 'bicubic');
            
            yfin(:,:,i) = yzdif(:,:,i) + imresize(ze_up(:,:,i), size(y(:,:,i)), 'bicubic'); % Combine denoised lower res image with denoised original image
        end
        PSNR_sub = -10*log10(var( (x(:)-yfin(:))));
        
        if PSNR_sub > psnrOut
            ze = yfin;
            PSNRbest(iter) = PSNR_sub;
        else
            PSNRbest(iter) = psnrOut;
        end
    else
        [sigma_hat] = noise_estimation(y, mode, sigma, 4, 1); % Estimate noise level; mode = 1 uses the true noise variance
        
        [ze, psnrOut] = denfun_global(x, y, sigma_hat, pdref, ps, Kref, mixture, dis, class, supportvar); % Denoise input image
                
        PSNR_sub = -10*log10(var( (x(:)-ze(:))));
    end
    
    psnrOut = -10*log10(var( (x(:) - ze(:))));
    PSNRbest(iter) = psnrOut;    
    
    if iter ~= 1
        
        final = 1/(1-ro)*ze - (ro/(1-ro))*prevze;
        PSNR_estimatefinal = -10*log10(var( (x(:) - final(:))));
        
        if PSNR_estimatefinal - psnrOut > 1e-2
            ze = final;
            PSNRbest(iter) = PSNR_estimatefinal;
        end
        cont = 0;
        
    end
    
    [mssim, ssim_map] = ssim(x, ze, [0.01 0.03], fspecial('gaussian', 11, 1.5), 1);
    
    prevze = ze;
        
    fprintf('Noisy PSNR = %0.5g; MMSE PSNR = %0.5g; SSIM = %0.5g\n', ...
    PSNR_noisy, PSNRbest(iter), mssim);
    if iter > 2
        cont = 0;
        break
    end
    
    iter = iter + 1;
end





