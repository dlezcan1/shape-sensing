function pos = position_kcp_3layer(B, Binv, kc_init1, kc_init2, kc_init3, w_init, L_init, ...
    L, s_crit1, s_crit2, p1, p2, p3)

%% Instantiations
    if s_crit1 >= L % 1-layer insertion
        pos = position_kcp(B, Binv, kc_init1, L_init, L, p1);
        
    elseif s_crit1 < L && L <= s_crit2 % 2-layer insertion
        pos = position_kcp_2layer(B, Binv, kc_init1, kc_init2, w_init, L_init, L, s_crit1, p1, p2);
        
    else % 3-layer insertion

        kc1 = kappa_c_p(kc_init1, L_init, s_crit1, p1);
        kc2 = kappa_c_p(kc_init2, L_init, s_crit2, p2);
        kc3 = kappa_c_p(kc_init3, L_init, L - s_crit2, p3);

        s = 0:.5:L;
    
%% Get the position
    wv = fn_intgEP_3layers_s_crit_Dimitri(w_init, kc1, kc2, kc3 ,...
        0, 0, 0.5, s_crit1, s_crit2, L, B, Binv);
    
    pos = wv2r(wv, L);   
    
    end
    
end
    