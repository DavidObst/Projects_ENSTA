n=10;// Nombre de pages

function show_adj(Adj, diameters)
  [lhs,rhs]=argn(0); 
  if rhs < 2 then diameters=30*ones(1,n);end
  graph=mat_2_graph(sparse(Adj),1,'node-node'); 
  graph('node_x')=300*cos(2*%pi*(1:n)/(n+1));
  graph('node_y')=300*sin(2*%pi*(1:n)/(n+1)); 
  graph('node_name')=string([1:n]);
  graph('node_diam')=diameters;
  //graph('node_color')= 1:n; 
  //show_graph(graph);
  rep=[1,1,1,1,2,2,2,2,2,2,2,2,2];
  plot_graph(graph,rep);
endfunction 


// Construction de la matrice de transition P 
// associée à une matrice d'adjacence.
// Pss: transition d'origine,
// P: matrice de google 
// z: vecteur de teleportation
// d: vecteur vaut 1 si le degré vaut zero et 0 sinon 

function [P, Pss, Pprim, d, z, alpha]=google(Adj)
  // <A completer>
  n=size(Adj,1);
  alpha=0.85;
  z=(1/n)*ones(1,n);
  d=zeros(1,n);
  for i=1:n
      if sum(Adj(i,:))==0
          d(i)=1;
      end
  end
  Pss=Adj;
  for i=1:n
      if (sum(Adj(i,:))>0)
          Pss(i,:)=Pss(i,:)/(sum(Adj(i,:)))
      end
  end
  Pprim=Pss;
  for i=1:n
      if (sum(Pss(i,:))==0)
          Pprim(i,:)=z;
      end
  end
  P=Pprim;
  P=alpha*P+(1-alpha)*((1/n)*ones(n,n));
endfunction
 
 
 
Adj=grand(n,n,'bin',1,0.2);//show_adj(Adj);
 
[P,Pss,Pprim,d,z,alpha]=google(Adj);
// verification que P est stochastique

sum(P,'c')


//le calcul de P′∗ x peut se faire en utilisant Pss′, d et z
// P=alpha*P1 + (1-alpha)* ones(n,1)*z
//P1=Pss+d'z
 x=rand(n,1)
y1=P'*x;
y2=alpha*Pss'*x+alpha*(d'*z)'*x+(1-alpha)*ones(n,1)*z*x;// A completer ...
y1-y2

//... = spec(...)
[eigenvec,eigenval] = spec(P');

// Recover the eigenvector corresponding to the eigenvalue 1
indice = find(abs(eigenval - 1) < 1e-12);
// A theorem ensures that a law pi exists
// However the function spec only yields an eigenvector proportional to pi
// Therefore we need to  normalize the obtained eigenvector corresponding to 1
// The normalization constant being the sum of the elements of the eigenvector
pi=eigenvec(:,indice)/sum(eigenvec(:,indice));
pi = pi';  // because pi is a LINE vector

xbasc();show_adj(Adj,int(300*pi)); // New plot of pagerank

// Question 5

function pi=pi_iterative()
    // We use value interation to calculate pi in another fashion
    p=ones(n,1);
    while %t do
      pn = P'*p ;
      if norm(pn-p,%inf) < 10*%eps then break;end
      p = pn ;
    end
    pi=(p/sum(p))';
endfunction
  
pi2=pi_iterative();
clean(pi*P-pi);

//Verification that the two invariants are indeed the same
abs(pi-pi2) < 1e-12

// Question 6 //
function pi=pi_iterative_sparse()
    p=ones(n,1);
    while %t do
      pn = alpha*Pss'*p+alpha*(d'*z)'*p+(1-alpha)*ones(n,1)*z*p ;
      if norm(pn-p,%inf) < 10*%eps then break;end
      p = pn ;
    end
    pi=(p/sum(p))';
endfunction
  
pi=pi_iterative_sparse();
clean(pi*P-pi)

// Verification 
abs(pi-pi2) < 1e-12


// maximization du << page rank >>
m=n/2;
p = 2;

// les liens des m premiers noeuds peuvent etre changés 
// et les m-premiers liens sont interdits i.e dans la 
// matrice d'adjacence on change Adj(1:m,m+1:$) 
// get_controls renvoit une matrix ncontxn 
// chaque ligne représentant un choix possible pour n-liens 

function controls=get_control(n)
  if n==1 then controls=[0;1];
  else
    A=get_control(n-1);
    controls = [zeros(size(A,'r'),1),A;ones(size(A,'r'),1),A];
  end
endfunction

controls = get_control(m);
// 2^5 controles 
ncont = size(controls,'r');

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pi=pi_iterative_sparse0(Pss,z,d,z)
    p=ones(n,1);
    while %t do
      pn = alpha*Pss'*p+alpha*(d'*z)'*p+(1-alpha)*ones(n,1)*z*p ;
      if norm(pn-p,%inf) < 10*%eps then break;end
      p = pn ;
    end
    pi=(p/sum(p))';
endfunction

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pi_opt,cost_opt,u_opt] = ranking(controls,p,m,Adj)
    Q = size(controls,'r');
    
    cost_opt = -10^5; //Optimal cost
    pi_opt = []; //Optimal pagerank
    u_opt = zeros(p,m); //Optimal control ;
    
    // Control other the pages m+1:n
    // Thus we can only change this submatrix
    // The adjacent matrix is initialized randomly
    
    Adj_ij = Adj;
    for i=1:Q
        for j=1:Q
            // We have access only to the pages 1->p (p=2 here)
            // The  constraints are that we cannot modify the m first pages.
            // However we can modify the m+1 -> n others
            // The modifications are the links to other pages which are shown or not
            // Those correspond to all the possible combinations 'controls' generated previously.
            Adj_ij(1,m+1:$) = controls(i,:);
            Adj_ij(2,m+1:$) = controls(j,:);
            
            // Generation of the P matrix as previously //
            [P,Pss,Pprim,d,z,alpha]=google(Adj_ij);
            
            // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
            // Calculation of the google ranking via the previous functions // 
            pi_ij = pi_iterative_sparse0(Pss,z,d,z);
            
            // Calculation of the cost : only the first m pages
            cost_ij = sum(pi_ij(1:m));
            
            if cost_ij>cost_opt
               cost_opt = cost_ij;
               pi_opt = pi_ij;
               u_opt(1,:) = controls(i,:); 
               u_opt(2,:) = controls(j,:);
            end
        end    
    end
endfunction
