function s_crit = determine_s_crit_1layer(kc1, w_init, z_crit, L, s0, ds, B, Binv)
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
    w_init = [kc1; 0; 0];
end
w1 = fn_intgEP_1layer_Dimitri(kc1,w_init,0,L,s0,ds,B,Binv);
p1 = wv2r(w1, L);
ix_crit_v = find(abs(p1(3,:) - z_crit) <= ds/2); 
if isempty(ix_crit_v)
    s_crit = -1;
    
else
    ix_crit = ix_crit_v(1);
    s_crit = s(ix_crit);
    
end

end