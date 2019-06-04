function [u] = upsamp(im)

[h, w] = size(im);

u1 = im(1:h/2, 1:w/2);
u2 = fliplr(im(1:h/2, w/2 + 1:end));
u3 = flipud(im(h/2 + 1:end, 1:w/2));
u4 = rot90(im(h/2 + 1:end, w/2 + 1:end),2);



for i = 0:h-1
    for j = 0:w-1
        p1 = floor((i)/2)+1;
        p2 = floor((i-1)/2)+1;
        q1 = floor((j)/2)+1;
        q2 = floor((j-1)/2)+1;
        
        if i == 0 && j == 0
            u(i+1,j+1) = u1(p1, q1);
        elseif i == 0
            u(i+1,j+1) = 0.5*(u1(p1, q1) + u2(p1, q2));
        elseif j == 0
            u(i+1,j+1) = 0.5*(u1(p1, q1) + u3(p2, q1));
        else
            u(i+1,j+1) = 0.25*(u1(p1, q1) + u2(p1, q2) + u3(p2, q1) + u4(p2, q2));
        end
    end
end

%figure, imagesc(u), colormap gray
