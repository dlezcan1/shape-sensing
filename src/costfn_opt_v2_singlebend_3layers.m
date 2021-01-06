function y = costfn_opt_v2_singlebend_3layers(eta,kc_gel,data,s_m,s0,ds,theta0,z1_crit,z2_crit,N,B,Binv,scalef)
%
% layer: gel-meat-gel
% kc for gel is obtained from the previous simulation.
% kc for meat is included in the optimization variables.
%
% - written by Jin Seob (Jesse) Kim

N_meas = size(data,2);
s_index_meas = s_m;

w_init = eta(1:3);
kc1 = kc_gel;
kc2 = eta(4);
kc3 = kc_gel;

[wv,pmat,Rmat] = fn_intgEP_v1_3layers(w_init,kc1,kc2,kc3,z1_crit,z2_crit,theta0,s0,ds,N,B,Binv);

% % include zero torsion
% yv = wv(1:3,s_index_meas) - data;
% y = norm(yv,'fro')^2*scalef;

% exclude torsion
yv = wv(1:2,s_index_meas) - data(1:2,:);
y = norm(yv,'fro')^2*scalef;


% L = (N-1)*ds; % in mm
% s = [0:ds:L];
% 
% % configuration
% wv = zeros(3,N);
% wv(:,1) = w_init;
% 
% Rmat = zeros(3,3,N);
% Rmat(:,:,1) = Rot_x(theta0);%eye(3);
% 
% pmat = zeros(3,N);
% for i = 2:N
%     z = pmat(3,i-1);
%     if z < z_crit
%         k0 = kc1*(1 - s(i-1)/L)^2;
%         w0 = [k0;0;0];
% 
%         k0prime = -2*kc1/L*(1 - s(i-1)/L);
%         w0prime = [k0prime;0;0];
%     else
%         k0 = kc2*(1 - s(i-1)/L)^2;
%         w0 = [k0;0;0];
% 
%         k0prime = -2*kc2/L*(1 - s(i-1)/L);
%         w0prime = [k0prime;0;0];
%     end
%     
%     if i == 2
%         wv(:,2) = w_init + ds*(w0prime - Binv*cross(w_init,B*(w_init - w0)));
%     else
%         wv(:,i) = wv(:,i-2) + 2*ds*(w0prime - Binv*cross(wv(:,i-1),B*(wv(:,i-1) - w0)));
%     end
%     
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
% 
% wv = fn_intgEP_v1(w_init,w0,w0prime,0,ds,N,B,Binv);
% 
% % % include zero torsion
% % yv = wv(1:3,s_index_meas) - data;
% % y = norm(yv,'fro')^2*scalef;
% 
% % exclude torsion
% yv = wv(1:2,s_index_meas) - data(1:2,:);
% y = norm(yv,'fro')^2*scalef;


% y1 = 0;
% for i = 1:3
%     y1 = y1 + norm(yv(:,i),2)^2;
% end

