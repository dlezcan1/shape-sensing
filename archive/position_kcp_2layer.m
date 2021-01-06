function pos = position_kcp_2layer(B, Binv, kc_init1, kc_init2, w_init, L_init, L, s_crit, p1, p2)
%% Instantiations
    if s_crit >= L || s_crit < 0
        pos = position_kcp(B, Binv, kc_init1, w_init, L_init, L, 0.492);
        
    elseif s_crit == 0
        pos = position_kcp(B, Binv, kc_init2, w_init, L_init, L, 0.492);
                
    else
%         kc1 = kappa_c_p(kc_init1, L_init, s_crit, p1);
%         kc2 = kappa_c_p(kc_init2, L_init, L-s_crit, p2);
        kc1 = kappa_c_p(kc_init1, L_init, L, p1);
        kc2 = kappa_c_p(kc_init2, L_init, L, p2);
    
    
    s = 0:.5:L;
    
    if L ~= L_init && length(w_init) == 3
        % w_init scaling
        r1 = kappa_c_p(1, L_init, s_crit, p1);
        r2 = kappa_c_p(1, L_init, L-s_crit, p2);   
%         r = (kc1*r1 + kc2*r2)/(kc1+kc2); % kc weighted average
%         r = (s_crit*r1 + (L-s_crit)*r2)/(L); % l weighted avg.
%         r = (kc1*s_crit*r1 + kc2*(L-s_crit)*r2)/(kc1*s_crit + kc2*(L-s_crit)); % kc-l weighted avg.
        
%         w_init = r * w_init;   
%         w_init = w_init * r1; % using p1 and 1st layer
%         w_init = w_init * r2; % using p2 and 2nd layer 
        w_init = w_init * kappa_c_p(1, L_init, L, p1); % using p1 and entire insertion (BEST)
        w_init = w_init * kappa_c_p(1, L_init, L, p2); %y using p2 and entire insertion
    
    end % if-else w_init
    
%% Get the position
    wv = fn_intgEP_2layers_s_crit_Dimitri(kc1, kc2, w_init, s_crit, L, 0, 0.5, B, Binv);
    
    pos = wv2r(wv, L);
    
    end % if-else insertino case
    
    