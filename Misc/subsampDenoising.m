function [ze_up] = subsampDenoising(y, x, sigma, K)

[~, ~, B] = size(y);

for i = 1:B
    ysub(:,:,i) = subsamp(y(:,:,i));
    xsub(:,:,i) = subsamp(x(:,:,i));
end

if sigma < 50
    pd = 6;
else
    pd = 8;
end

if B ~= 1
    pd = 4; 
end

sigma_hat = (sigma/255)/2;
[ze] = denfun_global(xsub, ysub, sigma_hat, pd, 1, K, '', 0, 0, []);

for i = 1:B
    ze_up(:,:,i) = upsamp(ze(:,:,i));
end

