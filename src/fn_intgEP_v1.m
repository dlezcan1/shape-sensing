function [wv] = fn_intgEP_v1(w_init,w0,w0prime,s0,ds,N,B,Binv)
%
% integration of Euler-Poincare equation for one layer
% (used in main_singlebend_v2_tissue_90total.m)
% w0: intrinsic curvature
% w0prime: d omega_0/ds
% theta0: initial rotation angle
% s0: initial arclength
% ds: delta s
% N: number of points in s
% B, Binv: stiffness matrix and its inverse
%
% - written by Jin Seob (Jesse) Kim

L = s0 + (N-1)*ds; % in mm
s = [s0:ds:L];

wv = zeros(3,N);
wv(:,1) = w_init;

% % linear
% k0 = kc*(1 - alpha*s);
% w0 = [k0;zeros(1,N);zeros(1,N)];
% 
% k0prime = -alpha*kc;
% w0prime = [k0prime*ones(1,N);zeros(1,N);zeros(1,N)];
% 
% % % exponential
% % k0 = kc*exp(-alpha*s);
% % w0 = [k0;zeros(1,N);zeros(1,N)];
% % 
% % k0prime = -alpha*k0;
% % w0prime = [k0prime;zeros(1,N);zeros(1,N)];

%% integration (trapezoidal method)
% i = 1
wv(:,2) = w_init + ds*(w0prime(:,1) - Binv*cross(w_init,B*(w_init - w0(:,1))));

% i > 1
for i = 2:N-1
    wv(:,i+1) = wv(:,i-1) + 2*ds*(w0prime(:,i) - Binv*cross(wv(:,i),B*(wv(:,i) - w0(:,i))));
end

% %% integration with 4th order Runge-Kutta method
% for i = 1:N-1
%     I1 = ds*fn_dwds(wv(:,i),w0(:,i),w0prime(:,i),B,Binv);
%     I2 = ds*fn_dwds(wv(:,i)+I1/2,w0(:,i),w0prime(:,i),B,Binv);
%     I3 = ds*fn_dwds(wv(:,i)+I2/2,w0(:,i),w0prime(:,i),B,Binv);
%     I4 = ds*fn_dwds(wv(:,i)+I3,w0(:,i),w0prime(:,i),B,Binv);
%     wv(:,i+1) = wv(:,i) + (I1 + 2*I2 + 2*I3 + I4)/6;
% end
%     
% %=================================
% % functions for the 4th order Runge-Kutta method
% function dwds = fn_dwds(w,w0,w0p,B,Binv)
% 
% dwds = w0p - Binv*cross(w,B*(w - w0));
    
    
    


