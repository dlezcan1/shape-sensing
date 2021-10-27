function [wv, pmat, Rmat] = fn_intgEP_1layer_Dimitri(kc,w_init, theta0, L,s0,ds,B,Binv)
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

    s = s0:ds:L;
    N = length(s);

    % intrinsic curvature
    k0 = kc*(1 - s/L).^2;
    w0 = [k0;zeros(1,N);zeros(1,N)];

    k0prime = -2*kc/L*(1 - s/L);
    w0prime = [k0prime;zeros(1,N);zeros(1,N)];

    % integrate
    if nargout == 1
        wv = fn_intgEP_w0_Dimitri(w_init, w0, w0prime, theta0, s0, ds, N, B, Binv);
    else
        [wv, pmat, Rmat] = fn_intgEP_w0_Dimitri(w_init, w0, w0prime, theta0, s0, ds, N, B, Binv);
    end
end