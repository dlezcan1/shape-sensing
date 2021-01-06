function [wv,fv] = fn_inextrod_w0(kc,p,f,m_ext,f_ext,B,Binv,ds,N)
% compute the needle trajectory given inputs
% this is to compute \omega_0 vector
% 
% - written by Jin Seob (Jesse) Kim

e3 = [0;0;1];

% f_ext = zeros(3,N);
% f_ext(2,:) = -p;
% f_ext(3,:) = -f;

wv = zeros(3,N);
wv(:,1) = [kc;0;0];

fv = zeros(3,N);
fv(:,1) = [0;p;f];%-f_ext(:,1);

for i = 2:N
%     % Euler integration
%     if i == 2
%         dfv = -cross(wv(:,i-1),fv(:,i-1)) - f_ext(:,i-1);
%         fv(:,i) = fv(:,i) + ds*dfv;
%         
%         wv(:,i) = wv(:,i-1) - ds*(Binv*cross(wv(:,i-1),B*wv(:,i-1)) + cross(e3,fv(:,i-1)));
%     else
%         dfv = -cross(wv(:,i-1),fv(:,i-1)) - f_ext(:,i-1);
%         fv(:,i) = fv(:,i-2) + 2*ds*dfv;
%         
%         wv(:,i) = wv(:,i-2) - 2*ds*(Binv*cross(wv(:,i-1),B*wv(:,i-1)) + cross(e3,fv(:,i-1)));
%     end
    
    % 4th order Runge-Kutta
    I1 = ds*fn_inextrodeq([wv(:,i-1);fv(:,i-1)],B,Binv,m_ext(:,i-1),f_ext(:,i-1));
    I2 = ds*fn_inextrodeq([wv(:,i-1);fv(:,i-1)]+I1/2,B,Binv,m_ext(:,i-1),f_ext(:,i-1));
    I3 = ds*fn_inextrodeq([wv(:,i-1);fv(:,i-1)]+I2/2,B,Binv,m_ext(:,i-1),f_ext(:,i-1));
    I4 = ds*fn_inextrodeq([wv(:,i-1);fv(:,i-1)]+I3,B,Binv,m_ext(:,i-1),f_ext(:,i-1));
    xv_next = [wv(:,i-1);fv(:,i-1)] + (I1 + 2*I2 + 2*I3 + I4)/6;
    wv(:,i) = xv_next(1:3);
    fv(:,i) = xv_next(4:6);
    
%     M1 = ds*fn_moment(wv(:,i-1),fv(:,i-1),B,Binv);
%     M2 = ds*fn_moment(wv(:,i-1)+M1/2,fv(:,i-1),B,Binv);
%     M3 = ds*fn_moment(wv(:,i-1)+M2/2,fv(:,i-1),B,Binv);
%     M4 = ds*fn_moment(wv(:,i-1)+M3,fv(:,i-1),B,Binv);
%     wv(:,i) = wv(:,i-1) + (M1 + 2*M2 + 2*M3 + M4)/6;
%     
%     F1 = ds*fn_force(wv(:,i-1),fv(:,i-1),f_ext(:,i-1));
%     F2 = ds*fn_force(wv(:,i-1),fv(:,i-1)+F1/2,f_ext(:,i-1));
%     F3 = ds*fn_force(wv(:,i-1),fv(:,i-1)+F2/2,f_ext(:,i-1));
%     F4 = ds*fn_force(wv(:,i-1),fv(:,i-1)+F3,f_ext(:,i-1));
%     fv(:,i) = fv(:,i-1) + (F1 + 2*F2 + 2*F3 + F4)/6;
end
    
%=================================
% functions for the 4th order Runge-Kutta method
function dxds = fn_inextrodeq(x,B,Binv,m_e,f_e)
e3 = [0;0;1];
w = x(1:3);
f = x(4:6);

dwds = -Binv*(cross(w,B*w) + cross(e3,f) - m_e);
dfds = -cross(w,f) - f_e;

dxds = [dwds;dfds];

% function dwds =fn_moment(w,f,B,Binv)
% e3 = [0;0;1];
% dwds = - Binv*(cross(w,B*w) + cross(e3,f));
% 
% function dfds = fn_force(w,f,f_e)
% dfds = - cross(w,f) - f_e;