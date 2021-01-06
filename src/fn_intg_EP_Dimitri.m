function [wv] = fn_intg_EP_Dimitri(w_init, w0, w0prime, L, s, ds, B, Binv)
%
% integration of Euler-Poincare equation for one layer
% for variable lengths
%
% w0: intrinsic curvature
% w0prime: d omega_0/ds
% theta0: initial rotation angle
% L: the length of the needle insertd
% s0: initial arclength
% ds: delta s
% N: number of points in s
% B, Binv: stiffness matrix and its inverse
%
% - written by Jin Seob (Jesse) Kim
% - edited by Dimitri Lezcano

N = length(s);
wv = zeros(3,N);
wv(:,1) = w_init;

%% integration (trapezoidal method)
% i = 1
wv(:,2) = w_init + ds*(w0prime(:,1) - Binv*cross(w_init,B*(w_init - w0(:,1))));

% i > 1
for i = 2:N-1
    wv(:,i+1) = wv(:,i-1) + 2*ds*(w0prime(:,i) - Binv*cross(wv(:,i),B*(wv(:,i) - w0(:,i))));
end

    
    


