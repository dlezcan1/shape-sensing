function wv = fn_intgEP_2layers_s_crit_Dimitri(kc1, kc2, w_init, s_crit, L, s0, ds, B, Binv)
%
% integration of Euler-Poincare equation for two layers
% kc1, kc2: kappa_c values for the layers (assumed to have kc-p implemented
% s_crit: arclength of the boundary between two layers (fixed)
% s0: initial arclength
% ds: delta s
% L: Length of insertion
% B, Binv: stiffness matrix and its inverse
%
% returns wv and s_crit
%
% - written by Jin Seob (Jesse) Kim
% - edited by Dimitri Lezcano

s = [s0:ds:L]; % mm
N = length(s);

if s_crit >= L
    wv = fn_intgEP_1layer_Dimitri(kc1, w_init, 0, L, 0, ds, B, Binv);

else
    %% find s_crit
    [s_crit_found,ix_crit] = min(abs(s - s_crit));

    s1 = s(1:ix_crit);
    s2 = s(ix_crit+1:N);

    %% EP integration parameters
    % kc1, kc2 input on it's own

    % % intrinsic curvature
    % % kc's measured by 90 mm insertion, so proportional consideration
    % kc1 = kc1*s_crit/L; % linear dependence of kappa_c 
    % kc2 = kc2*(L-s_crit)/L; % linear dependence of kappa_c

    k0 = zeros(1,N);
    k0(1:ix_crit) = kc1*(s_crit - s1).^2/L^2 + kc2*(1 - s_crit/L)*(1 + s_crit/L - 2*s1/L);
    k0(ix_crit+1:N) = kc2*(1 - s2/L).^2;
    w0 = [k0;zeros(1,N);zeros(1,N)];

    if length(w_init) ~= 3
        w_init = w0(:,1);
    end

    k0prime = zeros(1,N);
    k0prime(1:ix_crit) = -2*kc1/L^2*(s_crit - s1) - 2*kc2/L*(1 - s_crit/L);
    k0prime(ix_crit+1:N) = -2*kc2/L*(1 - s2/L);
    w0prime = [k0prime;zeros(1,N);zeros(1,N)];

    %% integrate to calculate omega vector
    wv = zeros(3,N);
    wv(:,1) = w_init;
    for i = 1:N-1
        if i == 1
            wv(:,i+1) = w_init + ds*(w0prime(:,1) - Binv*cross(w_init,B*(w_init - w0(:,1))));
        else
            wv(:,i+1) = wv(:,i-1) + 2*ds*(w0prime(:,i) - Binv*cross(wv(:,i),B*(wv(:,i) - w0(:,i))));
        end
    end
end

% %% orientation and position
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