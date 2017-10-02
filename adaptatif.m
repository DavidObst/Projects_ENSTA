function [xi_adapt,w]= adaptatif(delta,sigma_INS,sigma_BAR,sigma_ALT,r_0,v_0,r_INS,v_INS,T,N,map,h_ALT,c)
%% Initialisation%%

    xi_adapt = zeros(4,N,T);
    w = (1/N)*ones(N,T);
    
    for i=1:N
        xi_adapt(1:2,i,1) = r_0';
        xi_adapt(3:4,i,1) = v_0';
    end
    
    N_eff = 1/dot(w(:,1),w(:,1));
        
    for t=2:T
        
        if N_eff <= c*N
            xi_hat = fct_multi(xi_adapt(:,:,t-1),w(:,t-1)',N);

            %dr = xi_hat(1:2,:,t-1) - r_INS;
            %dv = xi_hat(3:4,:,t-1) - v_INS;

            for i=1:N
                dr = xi_hat(1:2,i) - r_INS(:,t-1);
                dv = xi_hat(3:4,i) - v_INS(:,t-1);

                dr = dr';
                dv = dv';

                %dx_temporaire = xi_sir(1:2,i,t)';
                %dv_temporaire = xi_sir(3:4,i,t)';

                [dx_temporaire,dv_temporaire] = h_update(dr,dv,sigma_INS,delta);

                xi_adapt(1:2,i,t) = r_INS(:,t) + dx_temporaire';
                xi_adapt(3:4,i,t) = v_INS(:,t) + dv_temporaire';

                %%%%%%%% MAJ DES POIDS %%%%%%%%

                [x,y] = coord(xi_adapt(1:2,i,t),map);
                w(i,t) = gaussien(h_ALT(t),map(x,y),sigma_BAR^2+sigma_ALT^2);
            end
            
        else
                xi_hat=xi_adapt(:,:,t-1);    
            for i=1:N

                dr = xi_hat(1:2,i) - r_INS(:,t-1);
                dv = xi_hat(3:4,i) - v_INS(:,t-1);

                dr = dr';
                dv = dv';

                [dx_temporaire,dv_temporaire] = h_update(dr,dv,sigma_INS,delta);

                xi_adapt(1:2,i,t) = r_INS(:,t) + dx_temporaire';
                xi_adapt(3:4,i,t) = v_INS(:,t) + dv_temporaire';

                %%%%%%%% MAJ DES POIDS %%%%%%%%

                [x,y] = coord(xi_adapt(1:2,i,t),map);
                w(i,t) = w(i,t-1)*gaussien(h_ALT(t),map(x,y),sigma_BAR^2+sigma_ALT^2);
            end
        end
        
        s = sum(w(:,t));
        w(:,t) = w(:,t)/s;
        
        N_eff = 1/dot(w(:,t),w(:,t));
    end
end