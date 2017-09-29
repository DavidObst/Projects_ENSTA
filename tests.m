N = 600;
[swag,swag2] = SIR(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT);

test_hat = zeros(4,T);

for t=1:T
   test_hat(:,t) = swag(:,:,t)*swag2(:,t); 
end

[swag,swag2] = SIR(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT);

test_hat2 = zeros(4,T);

for t=1:T
   test_hat2(:,t) = swag(:,:,t)*swag2(:,t); 
end

c = 100;
[swag,swag2] = adaptatif(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT,c);

test_hat3 = zeros(4,T);

for t=1:T
   test_hat3(:,t) = swag(:,:,t)*swag2(:,t); 
end

map = load('nct/mnt.data');
[N1 N2] = size(map);

imagesc([X1MIN X1MAX],[X2MIN X2MAX],transpose(map));
hold on;
colormap('default');
axis square;
axis off;

load('nct/traj.mat','rtrue','vtrue');
nmax = size(rtrue,2);

plot(rtrue(1,:),rtrue(2,:),'r-');
plot(r_INS(1,:),r_INS(2,:),'m-');
plot(test_hat(1,:),test_hat(2,:),'b');
plot(test_hat2(1,:),test_hat2(2,:),'g');
plot(test_hat3(1,:),test_hat3(2,:),'black');
legend('True','Non-filtered','SIR','SIS','Adapt')