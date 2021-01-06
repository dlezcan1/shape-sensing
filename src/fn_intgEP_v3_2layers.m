function [wv,pmat,Rmat] = fn_intgEP_v3_2layers(w_init,kc1,kc2,z_crit,theta0,s0,ds,N,B,Binv)
%
% integration of Euler-Poincare equation for two layers
% kc1, kc2: kappa_c values for the layers
% z_crit: z-coordinate of the boundary between two layers
% theta0: initial rotation angle
% s0: initial arclength
% ds: delta s
% N: number of points in s
% B, Binv: stiffness matrix and its inverse
%
% - written by Jin Seob (Jesse) Kim
% - edited by Dmitri Lezcano

L = s0 + (N-1)*ds; % in mm
s = [s0:ds:L];

% determine s_crit (90 mm)
[w1,p1] = fn_intgEP_v1_1layer(w_init,kc1,theta0,s0,ds,N,B,Binv);
ix_crit_v = find(abs(p1(3,:) - z_crit) <= ds/2); 
ix_crit = ix_crit_v(1);
s_crit = s(ix_crit);

s1 = s(1:ix_crit);
s2 = s(ix_crit+1:N);

% intrinsic curvature
k0 = zeros(1,N);
k0(1:ix_crit) = kc1*(s_crit - s1).^2/L^2 + kc2*(1 - s_crit/L)*(1 + s_crit/L - 2*s1/L);
k0(ix_crit+1:N) = kc2*(1 - s2/L).^2;
w0 = [k0;zeros(1,N);zeros(1,N)];

k0prime = zeros(1,N);
k0prime(1:ix_crit) = -2*kc1/L^2*(s_crit - s1) - 2*kc2/L*(1 - s_crit/L);
k0prime(ix_crit+1:N) = -2*kc2/L*(1 - s2/L);
w0prime = [k0prime;zeros(1,N);zeros(1,N)];

% integrate to calculate omega vector
wv = zeros(3,N);
wv(:,1) = w_init;
for i = 1:N-1
    if i == 1
        wv(:,i+1) = w_init + ds*(w0prime(:,1)-Binv*cross(w_init,B*(w_init - w0(:,1))));
    else
        wv(:,i+1) = wv(:,i-1) + 2*ds*(w0prime(:,i)-Binv*cross(wv(:,i),B*(wv(:,i)-w0(:,i))));
    end
end

% orientation and position
Rmat = zeros(3,3,N);
Rmat(:,:,1) = Rot_x(theta0);

pmat = zeros(3,N);
for i = 2:N
    % orientation
    W = matr(1/2*(wv(:,i-1) + wv(:,i)));
    Rmat(:,:,i) = Rmat(:,:,i-1)*expm(ds*W);

    % position
    e3vec = squeeze(Rmat(:,3,1:i));
    if i == 2
        pmat(:,i) = pmat(:,i-1) + squeeze(Rmat(:,3,i))*ds;
    else
        pmat(:,i) = Simpson_vec_int(e3vec,ds);
    end        
end
