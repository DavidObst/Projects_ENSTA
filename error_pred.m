function [ e ] = error_pred( rtrue,rhat,T )
%% average L2 error 
e=0;
for t=1:T
   e=e+ sqrt( dot(rtrue(:,t)-rhat(:,t),rtrue(:,t)-rhat(:,t)) );
end
e=e/T;
end

