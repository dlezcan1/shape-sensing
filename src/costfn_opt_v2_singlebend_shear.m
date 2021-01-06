function y = costfn_opt_v2_singlebend_shear(eta,data,s_m,ds,N,K,Kinv,scalef)

N_meas = size(data,2);
s_index_meas = s_m;

xi_init = eta(1:5);
k0 = eta(6);
v0y = eta(7);

L = (N-1)*ds; % in mm
s = [0:ds:L];

% % intrinsic deformation
% xi0 = [k0;0;0;0;v0y;1];

% integration of the E-P equation
xiv = fn_intgEP_v1_shear(xi_init,k0,v0y,0,ds,N,K,Kinv);

% % include zero torsion
% yv = wv(1:3,s_index_meas) - data;
% y = norm(yv,'fro')^2*scalef;

% exclude torsion
yv = xiv(1:2,s_index_meas) - data(1:2,:);
y = norm(yv,'fro')^2*scalef;


% y1 = 0;
% for i = 1:3
%     y1 = y1 + norm(yv(:,i),2)^2;
% end

