function [wv,pmat,Rmat] = fn_intgEP_w0_Dimitri(w_init, w0, w0prime,theta0,s0,ds,N,B,Binv)
%
% integration of Euler-Poincare equation for a homogeneous layer
% w0: the intrinsic curvature deformation vector ( open for flexibility )
% w0prime: the derfivative of w0
% theta0: initial rotation angle
% s0: initial arclength
% ds: delta s
% N: number of points in s
% B, Binv: stiffness matrix and its inverse
%
% - written by Jin Seob (Jesse) Kim
% - edited by Dimitri Lezcano

L = s0 + (N-1)*ds; % in mm
s = [s0:ds:L];

wv = zeros(3,N);
wv(:,1) = w_init;

% integrate to calculate omega vector
for i = 1:N-1
    if i == 1
        wv(:,i+1) = w_init + ds*(w0prime(:,1) - Binv*cross(w_init,B*(w_init - w0(:,1))));
    else
        wv(:,i+1) = wv(:,i-1) + 2*ds*(w0prime(:,i) - Binv*cross(wv(:,i),B*(wv(:,i) - w0(:,i))));
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
