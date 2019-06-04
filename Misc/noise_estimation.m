function [sigma_hat] = noise_estimation(image, mode, sigma, size, level)

if size > 2
    if mode == 0
        sigma_hat = (sigma/255)/1.4^(level-1);
    elseif mode == 1
        [sigma_hat] = NoiseLevel(image,size);
    elseif mode == 2
        hh = (image(1:length(image)-1,:)-image(2:length(image),:))'/sqrt(2);
        sigma_hat = mad(hh(:),1)/0.6745;
    elseif mode == 3
        hh = (image(1:length(image)-1,:)-image(2:length(image),:))'/sqrt(2);
        sigma_hat = mad(hh(:),0)*1.253;
    elseif mode == 4
        [sigma_hat] = estimateNoiseSDUsingKurts(image,size);
    else
        fprintf('Invalid parameter: mode\n')
    end
    if sigma_hat < 6e-3
        [sigma_hat] = noise_estimation( image , 1 , sigma , size-1 , 1);
    end
else
    %fprintf('Invalid parameter: size\n')
    [sigma_hat] = 6e-3;
end



% fprintf('Noise standard deviation: %d\n', sigma_hat*255)
