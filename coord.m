function [i,j] = coord(u,map)
    [N1 N2] = size(map);
    X1MAX = 10^4;
    X2MAX = 10^4;
    X1MIN = -10^4;
    X2MIN = -10^4;
    
    i = floor((N1/(2*X1MAX))*u(1)-((N1/(2*X1MAX))*X1MIN));
    j = floor((N2/(2*X2MAX))*u(2)-((N2/(2*X2MAX))*X2MIN));
end