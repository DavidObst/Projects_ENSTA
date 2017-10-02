function q = gaussien( x,m,v2 )
%% density of x according to the gaussian law N(m,v2)
q=(1/sqrt(2*pi*v2))*exp(-0.5*(1/v2)*dot(x-m,x-m));

end

