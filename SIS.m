function [xi_sis,w]= SIS(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT)
    %% Initialisation%%

    xi_sis = zeros(4,N,T);
    w = (1/N)*ones(N,T);
    
    for i=1:N
        xi_sis(1:2,i,1) = r_0';
        xi_sis(3:4,i,1) = v_0';
    end
        
    for t=2:T
                
        for i=1:N
            dr = xi_hat(1:2,i) - r_INS(:,t-1);
            dv = xi_hat(3:4,i) - v_INS(:,t-1);
            
            dr = dr';
            dv = dv';
                        
            [dx_temporaire,dv_temporaire] = h_update(dr,dv,sigma_INS,delta);
            
            xi_sis(1:2,i,t) = r_INS(:,t) + dx_temporaire';
            xi_sis(3:4,i,t) = v_INS(:,t) + dv_temporaire';
            
            %%%%%%%% MAJ DES POIDS %%%%%%%%
            
            [x,y] = coord(xi_sis(1:2,i,t),map);
            w(i,t) = w(i,t-1)*gaussien(h_ALT(t),map(x,y),sigma_BAR^2+sigma_ALT^2);
        end
        
        s = sum(w(:,t-1));
        w(:,t) = w(:,t)/s;
    end
end
