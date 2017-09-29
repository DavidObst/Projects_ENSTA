
%% Parameters for the display of the map %%
X1MIN = -10^4;
X1MAX = 10^4;
X2MIN = -10^4;
X2MAX = 10^4;

%% Other Parameters %%

r_0 = [-6000,-2000];
v_0 = [120,0];
sigma_r0 = 100;
sigma_v0 = 10 ;
sigma_INS = 7;
sigma_ALT = 10;
sigma_BAR = 20;
T = 100;
delta = 1;

%% Load the map %%
map = load('nct/mnt.data');
[N1 N2] = size(map);

imagesc([X1MIN X1MAX],[X2MIN X2MAX],transpose(map));
hold on;
colormap('default');
axis square;
axis off;

%% Load real trajectory r_k and plot it on the map%%
load('nct/traj.mat','rtrue','vtrue');
nmax = size(rtrue,2);

plot(rtrue(1,:),rtrue(2,:),'r-');

%% Load acceleration measures from inertial system %%
load('nct/ins.mat','a_INS');


%% Plot of the trajectory estimated using INS measures only %%
r_INS(:,1) = r_0; v_INS(:,1) = v_0; % initialisation
for n=2:nmax
    r_INS(:,n) = r_INS(:,n-1)+delta*v_INS(:,n-1);
    v_INS(:,n) = v_INS(:,n-1)+delta*a_INS(:,n-1);
end
plot(r_INS(1,:),r_INS(2,:),'m-');

%% Plot du relief
figure
load('nct/alt.mat');
plot(h_ALT,'+');
r_true_ij=rtrue;
r_true_ij(1,:)=floor((N1/(2*X1MAX))*rtrue(1,:)-((N1/(2*X1MAX))*X1MIN));
r_true_ij(2,:)=floor((N2/(2*X2MAX))*rtrue(2,:)-((N2/(2*X2MAX))*X2MIN));

plot(0:T,h_ALT,'+r',0:T,map(r_true_ij(1,:),r_true_ij(2,:)),'b')

%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load real trajectory r_k and plot it on the map%%
load('nct/traj.mat','rtrue','vtrue');
nmax = size(rtrue,2);

figure
plot(rtrue(1,:),rtrue(2,:),'r-');

N = 10
[swag,swag2] = SIR(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT)

test_hat = zeros(4,T);

for t=1:T
   test_hat(:,t) = swag(:,:,t)*swag2(:,t); 
end

plot(test_hat(1,:),test_hat(2,:));
