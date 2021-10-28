function [wv,pmat,Rmat,s_crit] = fn_intgEP_zcrit_2layers_Dimitri(w_init,kc1,kc2,z_crit,theta0,s0,ds,N,B,Binv)
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
    % - edited by: Dimitri Lezcano

    L = s0 + (N-1)*ds; % in mm
    s = s0:ds:L;

    %% determine s_crit (90 mm)
    [wv1,pmat1, Rmat1] = fn_intgEP_v1_1layer(w_init,kc1,theta0,s0,ds,N,B,Binv);
    if max(pmat1(3,:)) <= z_crit % single layer insertion
        wv = wv1; pmat = pmat1; Rmat = Rmat1; s_crit = -1;
        return; 
    end
    ix_crit_v = find(abs(pmat1(3,:) - z_crit) <= ds/2); 
    ix_crit = ix_crit_v(1);
    s_crit = s(ix_crit);

    s1 = s(1:ix_crit);
    s2 = s(ix_crit+1:N);

    %% intrinsic curvature
    % kappa_0 and kappa_0' calculations
    k0 = zeros(1,N);
    k0(1:ix_crit) = kc1*(s_crit - s1).^2/L^2 + kc2*(1 - s_crit/L)*(1 + s_crit/L - 2*s1/L);
    k0(ix_crit+1:N) = kc2*(1 - s2/L).^2;
    w0 = [k0;zeros(2,N)];

    k0prime = zeros(1,N);
    k0prime(1:ix_crit) = -2*kc1/L^2*(s_crit - s1) - 2*kc2/L*(1 - s_crit/L);
    k0prime(ix_crit+1:N) = -2*kc2/L*(1 - s2/L);
    w0prime = [k0prime;zeros(2,N)];

    % integrate
    if nargout == 1
        wv = fn_intgEP_w0_Dimitri(w_init, w0, w0prime, theta0, s0, ds, N, B, Binv);
    else
        [wv, pmat, Rmat] = fn_intgEP_w0_Dimitri(w_init, w0, w0prime, theta0, s0, ds, N, B, Binv);
    end     
end
