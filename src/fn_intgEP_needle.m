function [wv] = fn_intgEP_needle(w_init,w0,w0prime,ds,N,B,Binv)
%
% integration of Euler-Poincare equation
% this will be general enough to include any w0(s)
% w0: intrinsic curvature
% w0prime: d omega_0/ds
% ds: delta s
% N: number of points in s
% B, Binv: stiffness matrix and its inverse
%
% - written by Jin Seob (Jesse) Kim

L = (N-1)*ds; % in mm
s = [0:ds:L];

wv = zeros(3,N);
wv(:,1) = w_init;

% %% Euler integration (trapezoidal)
% % i = 1
% wv(:,2) = w_init + ds*(w0prime(:,1) - Binv*cross(w_init,B*(w_init - w0(:,1))));
% 
% % i > 1
% for i = 2:N-1
%     wv(:,i+1) = wv(:,i-1) + 2*ds*(w0prime(:,i) - Binv*cross(wv(:,i),B*(wv(:,i) - w0(:,i))));
% end

%% 4th order Runge-Kutta
for i = 2:N    
    I1 = ds*fn_EPeq(wv(:,i-1),w0(:,i-1),w0prime(:,i-1),B,Binv);
    I2 = ds*fn_EPeq(wv(:,i-1)+I1/2,w0(:,i-1),w0prime(:,i-1),B,Binv);
    I3 = ds*fn_EPeq(wv(:,i-1)+I2/2,w0(:,i-1),w0prime(:,i-1),B,Binv);
    I4 = ds*fn_EPeq(wv(:,i-1)+I3,w0(:,i-1),w0prime(:,i-1),B,Binv);
    wv(:,i) = wv(:,i-1) + (I1 + 2*I2 + 2*I3 + I4)/6;
end
    
%% =================================
% functions for the 4th order Runge-Kutta method
function dwds = fn_EPeq(w,w0,w0p,B,Binv)

dwds = w0p - Binv*cross(w,B*(w - w0));

