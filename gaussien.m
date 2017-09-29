function q = gaussien(x,m,v2)
    C = sqrt(2*pi*v2);
    q = (1/C)*exp(-0.5/v2*dot(x-m,x-m));
end