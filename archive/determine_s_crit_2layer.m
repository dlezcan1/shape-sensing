function [s_crit1,s_crit2] = determine_s_crit_2layer(kc_init1, kc_init2, w_init, z_crit1, z_crit2,...
    L_init, L, s0, ds, B, Binv, p1, p2)
%
% integration of Euler-Poincare equation for two layers
% kc1: kappa_c values for the first layer
% z_crit: z-coordinate of the boundary between two layers
% theta0: initial rotation angle
% s0: initial arclength
% ds: delta s
% L: the length of insertion
% B, Binv: stiffness matrix and its inverse
%
% returns wv and s_crit
%
% - edited by Dimitri Lezcano

% determine s_crit

s = s0:ds:L;
if w_init ~= 3
    w_init = [kc_init1; 0; 0];
end

%% determine s_crit1 (90 mm)
% [w1, throw] = fn_intgEP_1layer_Dimitri(kc1, w_init, 0, L_init, s0, ds, B, Binv);
% p_1layer = wv2r(w1);
p_1layer = position_kcp(B, Binv, kc_init1, L_init, L_init, p1);
ix_crit_v1 = find(abs(p_1layer(3,:) - z_crit1) <= ds/2); 

if isempty(ix_crit_v1)
    s_crit1 = -1;
    s_crit2 = -1;
    
else
    ix_crit1 = ix_crit_v1(1);
    s_crit1 = s(ix_crit1);

    %% determine s_crit2
    p_2layer = position_kcp_2layer(B, Binv, kc_init1, kc_init2, w_init, L_init, L, s_crit1, p1, p2);

    ix_crit_v2 = find(abs(p_2layer(3,:) - z_crit2) <= ds/2); 

    if isempty(ix_crit_v2)
        s_crit2 = -1;

    else
        ix_crit2 = ix_crit_v2(1);
        s_crit2 = s(ix_crit2);

    end
    
end

end