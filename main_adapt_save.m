clear

N = 300; %% Nombre de particules
Nb = 10;
lambda = 20;
%c = [500,200,2];
%c = [300,400,100];
c = [500,400,5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SEQUENCE = './seq1/';
START = 1;
% charge le nom des images de la sequence
filenames = dir([SEQUENCE '*.png']);
filenames = sort({filenames.name});
T = length(filenames);
% charge la premiere image dans 'im'
tt = START;
im = imread([SEQUENCE filenames{tt}]);
% affiche 'im'
figure;
set(gcf,'DoubleBuffer','on');
imagesc(im);

disp('Cliquer 4 points dans l''image pour definir la zone a suivre.');
% on recupere la zone a tracker
zone = zeros(2,4);
compteur=1;
while(compteur ~= 5)
    [x,y,button] = ginput(1);
    zone(1,compteur) = x;
    zone(2,compteur) = y;
    text(x,y,'X','Color','r');
    compteur = compteur+1;
end
newzone = zeros(2,4);
newzone(1,:) = sort(zone(1,:));
newzone(2,:) = sort(zone(2,:));
% definition de la zone a tracker
% x haut gauche, y haut gauche, largeur, hauteur
zoneAT = zeros(1,4);
zoneAT(1) = newzone(1,1);
zoneAT(2) = newzone(2,1);
zoneAT(3) = newzone(1,4)-newzone(1,1);
zoneAT(4) = newzone(2,4)-newzone(2,1);
% affichage du rectangle
rectangle('Position',zoneAT,'EdgeColor','r','LineWidth',3);
 
littleim = imcrop(im,zoneAT(1:4));
[tmp,Cmap] = rgb2ind(littleim,Nb,'nodither');
littleim = rgb2ind(littleim,Cmap,'nodither');
histoRef = imhist(littleim,Cmap);
histoRef = histoRef/sum(histoRef);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialization of C %%%
C = diag(c);

 %% Init des particules %%
 Id = eye([3,3]);
 sigma0 = 0.2;
 %pos0 = [zoneAT(1) + zoneAT(3)/2; zoneAT(2)-zoneAT(4)/2]; %% Prendre le centre du rectangle
 rect0 = zeros(4,2);
 rect0(1,:) = [zoneAT(1),zoneAT(2)];
 rect0(2,:) = [zoneAT(1)+zoneAT(3),zoneAT(2)];
 rect0(3,:) = [zoneAT(1),zoneAT(2)+zoneAT(4)];
 rect0(4,:) = [zoneAT(1)+zoneAT(3),zoneAT(2)+zoneAT(4)];
 pos0 = 0.5*(rect0(1,:) + rect0(4,:));
 pos0 = pos0';

 
 s0 = 100; %% Pourcentage initial
 
 %% Coordonees : [x,y,s] %%
 %X0 = normrnd([pos0,s0],sigma0,[3,N]);
 X0 = [pos0;s0];
 Xi_k = zeros(3,N,T);
 Xi_k(:,:,1) = mvnrnd(X0,sigma0*Id,N)';
 w = (1/N)*ones(N,T); %% Poids uniformes intialement
 zone_tt = zeros(T,N,4);
 zone_hat = zeros(T,4);

 pos_hat = zeros(T,3);
 rect_hat = zeros(T,1);
 
 for i=1:N
     zone_tt(1,i,:) = zoneAT;
 end
 
 
 c = 0.3 ;
 N_eff = 1/dot(w(:,1),w(:,1));
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  
 while tt<T %%Tant que le film n est pas termine
     tt = tt+1;
     im = imread([SEQUENCE filenames{tt}]);
     
     if N_eff<= c*N
	 
         %% Lecture de l'image

        Xi_tempo = fct_multi(Xi_k(:,:,tt-1),w(:,tt-1)',N);

        %% Update at step k %%

         for i=1:N

            U = mvnrnd([0,0,0],C);
            Xi_k(:,i,tt) = Xi_tempo(:,i) + U';

            zone_tt(tt,i,:) = transforme(zoneAT,X0',Xi_k(:,i,tt)');


            impart = imcrop(im,zone_tt(tt,i,:));
             impart = rgb2ind(impart,Cmap,'nodither');
             
             
               if ~isempty(impart)
                 histo = imhist(impart,Cmap);
                 histo = histo/sum(histo);
             end

             if isempty(impart) %% If particle out of screen : set its weight to 0
                 w(i,tt) = 0;
             else
                 w(i,tt) = vraisemblance(histoRef,histo,lambda);
             end


         end
         
     else
        
         for i=1:N

            U = mvnrnd([0,0,0],C);
            Xi_k(:,i,tt) = Xi_k(:,i,tt-1) + U';
            
             zone_tt(tt,i,:) = transforme(zoneAT,X0',Xi_k(:,i,tt)');


            impart = imcrop(im,zone_tt(tt,i,:));
             impart = rgb2ind(impart,Cmap,'nodither');
             
             if ~isempty(impart)
                 histo = imhist(impart,Cmap);
                 histo = histo/sum(histo);
             end

             if isempty(impart) %% If particle out of screen : set its weight to 0
                 w(i,tt) = 0;
             else
                 w(i,tt) = w(i,tt-1)*vraisemblance(histoRef,histo,lambda);
             end


         end
         
     end
     
     w(:,tt) = w(:,tt)/sum(w(:,tt));
     
     N_eff = 1/dot(w(:,tt),w(:,tt));
     
     pos_hat(tt,1) = Xi_k(1,:,tt)*w(:,tt) ;
     pos_hat(tt,2) = Xi_k(2,:,tt)*w(:,tt) ;
     pos_hat(tt,3) = Xi_k(3,:,tt)*w(:,tt) ;
     
     zone_hat(tt,:) = transforme(zoneAT,X0',pos_hat(tt,:));
     
     %figure
     clf;
     imagesc(im)
     title(['N_{eff} =' num2str(N_eff) '; cN =' num2str(c*N) ])
     hold on
    %scatter(pos_hat(tt,1),pos_hat(tt,2));
    scatter(Xi_k(1,:,tt),Xi_k(2,:,tt));
    
    rectangle('Position',zone_hat(tt,:),'EdgeColor','r','LineWidth',3);
    %for i=1:N
     %   rectangle('Position',zone_tt(tt,i,:),'EdgeColor','blue','LineWidth',1);
    %end
    drawnow;
    
    tt

 end