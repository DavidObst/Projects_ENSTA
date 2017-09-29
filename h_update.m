function  [drn,dvn] = h_update(dr,dv,v,Delta)
    w = normrnd(0,v,1,2);
    drn = dr+Delta*dv;
    dvn = dv - Delta*w;
end

