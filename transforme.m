function [ zoneAT1 ] = transforme( zoneAT0,X0,X1 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
zoneAT1=zeros(1,4);
zoneAT1(1:2) = X1(1:2)-X0(1:2)+zoneAT0(1:2);
R=X1(3)/X0(3);
v=(zoneAT1(1:2)-X1(1:2))*(1-R);
zoneAT1(1:2)=zoneAT1(1:2)-v;
zoneAT1(3:4)=zoneAT0(3:4)*R;
end

