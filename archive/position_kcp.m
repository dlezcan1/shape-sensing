function p = position_kcp(B,Binv,kc_init, w_init, Linit,  L, p)
%% Instantiations
    kc = kappa_c_p(kc_init,Linit,L,p);
    s = 0:.5:L;
    
    if length(w_init) < 3
        w_init = [kc; 0; 0];
    end
    
    k0 = kc*(1 - s/L).^2;
    k0_prime = -2*kc/L*(1 - s/L);
    w0 = [k0; zeros(2,length(s))];
    w0_prime = [k0_prime; zeros(2,length(s))];
    
%% Get the position
    wv = fn_intg_EP_Dimitri(w_init, w0, w0_prime, L, s, 0.5, B, Binv);
    
    p = wv2r(wv, L);
    
    
end
    
    