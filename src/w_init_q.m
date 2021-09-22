function w_init = w_init_q ( w_init_ref, L_ref, L , q )
    w_init = w_init_ref .* (L_ref / L)^q;
    
end