function wv = fn_intgEP_1layer_Dimitri(kc,w_init, theta0, L,s0,ds,B,Binv)
%
% integration of Euler-Poincare equation for a homogeneous layer
% kc: kappa_c value for the layer
% theta0: initial rotation angle

% s0: initial arclength
% L: the length for insertion
% ds: delta s
% B, Binv: stiffness matrix and its inverse
%
% - written by Jin Seob (Jesse) Kim
% - edited by Dimitri Lezcano

if length(w_init) ~= 3
    w_init = [kc; 0; 0]; % ideal insertion
    
end

s = [s0:ds:L];
N = length(s);

wv = zeros(3,N);
wv(:,1) = w_init;

% intrinsic curvature
k0 = kc*(1 - s/L).^2;
w0 = [k0;zeros(1,N);zeros(1,N)];

k0prime = -2*kc/L*(1 - s/L);
w0prime = [k0prime;zeros(1,N);zeros(1,N)];

% integrate to calculate omega vector
for i = 1:N-1
    if i == 1
        wv(:,i+1) = w_init + ds*(w0prime(:,1) - Binv*cross(w_init,B*(w_init - w0(:,1))));
    else
        wv(:,i+1) = wv(:,i-1) + 2*ds*(w0prime(:,i) - Binv*cross(wv(:,i),B*(wv(:,i) - w0(:,i))));
    end
end

% % orientation and position
% Rmat = zeros(3,3,N);
% Rmat(:,:,1) = Rot_x(theta0);
% 
% pmat = zeros(3,N);
% for i = 2:N
%     % orientation
%     W = matr(1/2*(wv(:,i-1) + wv(:,i)));
%     Rmat(:,:,i) = Rmat(:,:,i-1)*expm(ds*W);
% 
%     % position
%     e3vec = squeeze(Rmat(:,3,1:i));
%     if i == 2
%         pmat(:,i) = pmat(:,i-1) + squeeze(Rmat(:,3,i))*ds;
%     else
%         pmat(:,i) = Simpson_vec_int(e3vec,ds);
%     end        
% end

end