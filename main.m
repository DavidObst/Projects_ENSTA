clear

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
figure
plot(0:T,h_ALT,'+r',0:T,map(r_true_ij(1,:),r_true_ij(2,:)),'b')
%plot(0:T,map(r_true_ij(1,:),r_true_ij(2,:)),'b')

%% Testes SIR
N = 2000;
[swag,swag2] = SIR(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT);

test_hat = zeros(4,T);

for t=1:T
   test_hat(:,t) = swag(:,:,t)*swag2(:,t); 
end

imagesc([X1MIN X1MAX],[X2MIN X2MAX],transpose(map));
hold on;
colormap('default');
axis square;
axis off;

plot(rtrue(1,:),rtrue(2,:),'r-');
plot(test_hat(1,:),test_hat(2,:),'b');

%% Testes SIS
N = 2000;
[swag,swag2] = SIS(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT);

test_hat = zeros(4,T);

for t=1:T
   test_hat(:,t) = swag(:,:,t)*swag2(:,t); 
end


imagesc([X1MIN X1MAX],[X2MIN X2MAX],transpose(map));
hold on;
colormap('default');
axis square;
axis off;

plot(rtrue(1,:),rtrue(2,:),'r-');
plot(test_hat(1,:),test_hat(2,:),'b');
%% Test adaptatif
N = 2000;
c=.5;
[swag,swag2] = adaptatif(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT,c);

test_hat3 = zeros(4,T);

for t=1:T
   test_hat3(:,t) = swag(:,:,t)*swag2(:,t); 
end

imagesc([X1MIN X1MAX],[X2MIN X2MAX],transpose(map));
hold on;
colormap('default');
axis square;
axis off;

plot(rtrue(1,:),rtrue(2,:),'r-');
plot(test_hat3(1,:),test_hat3(2,:),'b');

%% error tests SIR

N_set = [100:20:500,550:50:1000,5000];
N_essaie=50;
e_SIR=zeros(length(N_set),N_essaie);
e_mean_SIR=zeros(1,length(N_set));
e_var_SIR=zeros(1,length(N_set));
seed=rng;
for i=1:length(N_set)
    N=N_set(i);
    rng(seed); %set back to the same random queue 
    for j=1:N_essaie
        [swag,swag2] = SIR(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT);

        test_hat = zeros(4,T);

        for t=1:T
           test_hat(:,t) = swag(:,:,t)*swag2(:,t); 
        end

        rhat=test_hat(1:2,:);
        e_SIR(i,j)=error_pred( rtrue,rhat,T );
    end
    e_mean_SIR(i)=sum(e_SIR(i,:))/N_essaie; 
    e_var_SIR(i)=dot(e_SIR(i,:)-e_mean_SIR(i),e_SIR(i,:)-e_mean_SIR(i))/N_essaie;
end

%dans nct, j'ai ajouté les matrices des experiences avec N_essaye=50,
load('nct/e_SIR.mat');
for i=1:l
    e_mean_SIR(i)=sum(e_SIR(i,:))/N_essaie; 
    e_var_SIR(i)=dot(e_SIR(i,:)-e_mean_SIR(i),e_SIR(i,:)-e_mean_SIR(i))/N_essaie;
end

%l=length(N_set)+1; %full N_set
l=length(N_set); %doesn't show the N=5000 point
plot(N_set(1:l-1),e_mean_SIR(1:l-1),'b',...
    N_set(1:l-1),(e_mean_SIR(1:l-1)-(1.96*sqrt(e_var_SIR(1:l-1)/N_essaie))),'--r',...
    N_set(1:l-1),e_mean_SIR(1:l-1)+(1.96*sqrt(e_var_SIR(1:l-1)/N_essaie)),'--r');

%% error tests SIS

N_set = [100:20:500,550:50:1000,5000];
N_essaie=50;
e_SIS=zeros(length(N_set),N_essaie);
e_mean_SIS=zeros(1,length(N_set));
e_var_SIS=zeros(1,length(N_set));
seed=rng;
for i=1:length(N_set)
    N=N_set(i);
    rng(seed); %set back to the same random queue 
    for j=1:N_essaie
        [swag,swag2] = SIS(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT);
        
        test_hat = zeros(4,T);

        for t=1:T
           test_hat(:,t) = swag(:,:,t)*swag2(:,t); 
        end

        rhat=test_hat(1:2,:);
        e_SIS(i,j)=error_pred( rtrue,rhat,T );
    end
    e_mean_SIS(i)=sum(e_SIS(i,:))/N_essaie;
    e_var_SIS(i)=dot(e_SIS(i,:)-e_mean_SIS(i),e_SIS(i,:)-e_mean_SIS(i))/N_essaie;
end

%dans nct, j'ai ajouté les matrices des experiences avec N_essaye=50,
load('nct/e_SIS.mat');
for i=1:l
    e_mean_SIS(i)=sum(e_SIS(i,:))/N_essaie; 
    e_var_SIS(i)=dot(e_SIS(i,:)-e_mean_SIS(i),e_SIS(i,:)-e_mean_SIS(i))/N_essaie;
end

%l=length(N_set)+1; %full N_set
l=length(N_set); %doesn't show the N=5000 point
plot(N_set(1:l-1),e_mean_SIS(1:l-1),'b',...
    N_set(1:l-1),(e_mean_SIS(1:l-1)-(1.96*sqrt(e_var_SIS(1:l-1)/N_essaie))),'--r',...
    N_set(1:l-1),e_mean_SIS(1:l-1)+(1.96*sqrt(e_var_SIS(1:l-1)/N_essaie)),'--r');

%% error tests adaptatif

N_set = 100:50:1000;
C_set=0:0.2:1;
N_essaie=50;
e_adap=zeros(length(N_set),length(C_set),N_essaie);
e_mean_adap=zeros(length(C_set),length(N_set));
e_var_adap=zeros(length(C_set),length(N_set));
seed=rng;
for n=1:length(N_set)
    N=N_set(n);
    for c=1:length(C_set) 
        C=C_set(c);
        rng(seed); %set back to the same random queue 
        for j=1:N_essaie
            [swag,swag2] = adaptatif(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT,C);

            test_hat = zeros(4,T);

            for t=1:T
               test_hat(:,t) = swag(:,:,t)*swag2(:,t); 
            end

            rhat=test_hat(1:2,:);
            e_adap(n,c,j)=error_pred( rtrue,rhat,T );
        end
        e_mean_adap(c,n)=sum(e_adap(n,c,:))/N_essaie;
        e_var_adap(c,n)=dot(e_adap(n,c,:)-e_mean_adap(c,n),e_adap(n,c,:)-e_mean_adap(c,n))/N_essaie;
    end
end

%dans nct, j'ai ajouté les matrices des experiences avec N_essaye=50,
%si tu veux relance avec un plus grands N_essaye sur les pc de l'ensta 3:) 
load('nct/e_adap.mat');
for n=1:length(N_set)
    for c=1:length(C_set) 
        e_mean_adap(c,n)=sum(e_adap(n,c,:))/N_essaie;
        e_var_adap(c,n)=dot(e_adap(n,c,:)-e_mean_adap(c,n),e_adap(n,c,:)-e_mean_adap(c,n))/N_essaie;
    end
end
%manque d'afficher e_mean_adap et e_var_adap ... essaye avec HeatMap si tu peux :), 
%j'ai pas sue bien utiliser le parametre colormap
