function PSNRdb = PSNR(x, y, maxval, borders)
    if ~exist('borders', 'var'), borders = 0; end
    if ~exist('maxval', 'var'), maxval = 255; end
    
    xx=borders+1:size(x,1)-borders;
    yy=borders+1:size(x,2)-borders;
            
    PSNRdb = zeros(1,size(x,3));
    for fr=1:size(x,3) 
        err = x(xx,yy,fr) - y(xx,yy,fr);
        PSNRdb(fr) = 10 * log10((maxval^2)/mean(mean(err.^2)));    
    end
end