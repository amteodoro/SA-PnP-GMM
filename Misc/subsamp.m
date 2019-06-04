function [um] = subsamp(im)
% Sub-sampling from noise-clinic

[h, w] = size(im);

u1 = zeros(floor(h/2 - 1), floor(w/2 - 1));
u2 = u1;
u3 = u1;
u4 = u1;

for i = 0:floor(h/2-1)
    
    for j = 0:floor(w/2-1)
        
        u1(i+1,j+1) = 0.25*(im(2*i+1, 2*j+1) + im(2*i+1 + 1, 2*j+1) +im(2*i+1, 2*j+1 + 1) + im(2*i+1 + 1, 2*j+1 + 1));
        
        if j == floor(w/2 - 1)
            
            u2(i+1,j+1) = 0.5*(im(2*i+1, 2*j+1 + 1) + im(2*i+1 + 1, 2*j+1 + 1));
            
        else
            
            u2(i+1,j+1) = 0.25*(im(2*i+1, 2*j+1 + 1) + im(2*i+1 + 1, 2*j+1 + 1) +im(2*i+1, 2*j+1 + 2) + im(2*i+1 + 1, 2*j+1 + 2));
            
        end
        
        if i == floor(h/2 - 1)
            
            u3(i+1,j+1) = 0.5*(im(2*i+1 + 1, 2*j+1) + im(2*i+1 + 1, 2*j+1 + 1));
            
        else
            
            u3(i+1,j+1) = 0.25*(im(2*i+1 + 1, 2*j+1) + im(2*i+1 + 2, 2*j+1) +im(2*i+1 + 1, 2*j+1 + 1) + im(2*i+1 + 2, 2*j+1 + 1));
            
        end

        if i == floor(h/2 - 1) && j == floor(w/2 - 1)
            
            u4(i+1,j+1) = im(2*i+1 + 1, 2*j+1 + 1);
            
        elseif i == floor(h/2 - 1)
            
            u4(i+1,j+1) = 0.5*(im(2*i+1 + 1, 2*j+1 + 1) + im(2*i+1 + 1, 2*j+1 + 2));
            
        elseif j == floor(w/2 - 1)
            
            u4(i+1,j+1) = 0.5*(im(2*i+1 + 1, 2*j+1 + 1) + im(2*i+1 + 2, 2*j+1 + 1));
            
        else
            
            u4(i+1,j+1) = 0.25*(im(2*i+1 + 1, 2*j+1 + 1) + im(2*i+1 + 2, 2*j+1 + 1) + im(2*i+1 + 1, 2*j+1 + 2) + im(2*i+1 + 2, 2*j+1 + 2));
        end
    end
end

%subs = {u1; u2; u3; u4};

um = [u1, fliplr(u2); flipud(u3), rot90(u4,2)];

%figure, imagesc(um), colormap gray


% N = 2;
% 
% [ws, hs] = size(u1);
% 
% wm = N*ws;
% hm = N*hs;
% 
% for p = 0:N-1
%     for q = 0:N-1
%         for i = 1:hs-1
%             for j = 1:ws-1
%                 
%                 n = p*N + q + 1;
%                 
%                 if mod(p, 2) == 0
%                     di = i;
%                 else
%                     di = hs - 1 - i;
%                 end
%                 
%                 if mod(q, 2) == 0
%                     dj = j;
%                 else
%                     dj = ws - 1 - j;
%                 end
%                 
%                 um(di + p*hs + i, dj + q*ws + j) = subs{n}(i,j);
%                 
%             end
%         end
%     end
% end
%       
% 
% figure, imagesc(um), colormap gray

