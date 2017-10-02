function [ i,j ] = coord( u, map )
%% linear tranformation 
X1MIN = -10^4;
X1MAX = 10^4;
X2MIN = -10^4;
X2MAX = 10^4;
[N1 N2] = size(map);

i=floor((N1/(2*X1MAX))*u(1)-((N1/(2*X1MAX))*X1MIN));
j=floor((N2/(2*X2MAX))*u(2)-((N2/(2*X2MAX))*X2MIN));

%% case where simulated trajectory u outbounds the map
if j<1 
    j=1;
elseif j>N2 
    j=N2;
end

if i<1 
    i=1;
elseif i>N1 
    i=N1;
end

end

