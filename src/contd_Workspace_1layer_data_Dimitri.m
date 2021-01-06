%%
%
% continued from Workspace_1layer_data_Dimitri.m
% to test whether all are good.

% first without rotation
% p is given and q is computed 

% first needle trajectory at L = 90mm
%[p_ref, wv_ref, ~] = fn_intgEP_v1_1layer(w_init_ref,kc_s,theta0,0,ds,N_arclengths(1),B,Binv);
[p_ref, wv_ref] = integrate_twist(kc_s, L, N_arclengths(1), q, 0);

% needle trajectory without rotation
[p_org, wv_org] = integrate_twist(kc_s, L_pred, N_pred, q, 0); 

figure;
plot(squeeze(p_ref(3,:)), squeeze(p_ref(2,:)),'k-','LineWidth',2)
axis equal
axis([0 90 -15 15])
grid on

figure;
plot(squeeze(p_ref(3,:)), squeeze(p_ref(2,:)),'k-','LineWidth',2)
hold on
plot(squeeze(p_org(3,:)), squeeze(p_org(2,:)),'b-','LineWidth',2)
axis equal
axis([0 180 -40 40])
grid on

% rotate the needle at L = 90mm by 180 degrees
wv_rot = zeros(size(wv_org));
wv_rot(:,1:N_arclengths(1)) = wv_org(:,1:N_arclengths(1));
wv_rot(:,N_arclengths(1)+1:end) = Rot_z(pi)*wv_org(:,N_arclengths(1)+1:end);

p_rot = intg_MatDiffEq(wv_rot);

figure;
plot(squeeze(p_org(3,:)), squeeze(p_org(2,:)),'k-','LineWidth',2)
hold on
plot(squeeze(p_rot(3,:)), squeeze(p_rot(2,:)),'b-','LineWidth',2)
axis equal
axis([0 180 -40 40])
grid on

figure;
plot([0:ds:L_pred], wv_org(1,:))
hold on
plot([0:ds:L_pred], wv_rot(1,:))
axis([0 180 -0.005 0.005])
grid on














% function for integrating using a rotation at reference L
function [pmat, wv] = integrate_twist(kc, Li, Ni, q, theta_rot)
    global w_init_ref N_arclengths L p ds B Binv theta0
    
    % parameter scaling
    kc_i = kc*(L/Li)^p; % scale the kappa_c value\
    w_init = w_init_ref * (L/Li)^q;
    
    % arclength
    s = 0:ds:Li;
    
    % intrinsic curvature
    k0 = kc_i*(1 - s/Li).^2;
    w0 = [k0;zeros(1,Ni);zeros(1,Ni)];

    k0prime = -2*kc_i/Li*(1 - s/Li);
    w0prime = [k0prime;zeros(1,Ni);zeros(1,Ni)];

    % Rotate using rotation z matriz 
    Rz = Rot_z(theta_rot);
    w0(:,N_arclengths(1)+1:end) = Rz*w0(:,N_arclengths(1)+1:end);
    w0prime(:,N_arclengths(1)+1:end) = Rz*w0prime(:,N_arclengths(1)+1:end);

    [wv,pmat,~] = fn_intgEP_w0_Dimitri(w_init, w0, w0prime,theta0,0,...
        ds,Ni,B,Binv);
end


function [pmat, Rmat] = intg_MatDiffEq(wv)

global w_init_ref N_arclengths L p ds B Binv theta0

N = size(wv,2);
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

end