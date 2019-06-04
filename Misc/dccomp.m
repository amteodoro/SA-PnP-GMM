function [x_hat, patches] = dccomp(ac, dc, weights, pd, ps, d, dd)

global_dc = 0;

patches=bsxfun(@plus,ac,dc);

patches = min(max( patches , - global_dc) , 1 - global_dc);

x_hat_raw = col2imstep( patches.*(weights) , dd , [pd pd] , [ps ps]);
normalize = col2imstep(weights , dd , [pd pd] , [ps ps]);
x_hat1 = x_hat_raw ./ normalize;

x_hat = x_hat1( pd : ( d(1) + pd - 1) , pd : ( d(2) + pd - 1 ) ) + global_dc;

x_hat = min( max( x_hat , 0 ) , 1);