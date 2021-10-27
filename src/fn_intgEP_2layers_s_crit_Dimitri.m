function [wv, pmat, Rmat] = fn_intgEP_2layers_s_crit_Dimitri(kc1, kc2, w_init, s_crit, L, s0, ds, B, Binv)
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

    s = s0:ds:L; % mm
    N = length(s);

    if s_crit >= L % single-layer
        warning('Single-layer used in double-layer code.');
        if nargout == 1
            wv = fn_intgEP_1layer_Dimitri(kc1, w_init, 0, L, s0, ds, B, Binv);
        else
            [wv, pmat, Rmat] = fn_intgEP_1layer_Dimitri(kc1, w_init, 0, L, s0, ds, B, Binv);
        end

    else % double-layer
        %% find s_crit
        [~,ix_crit] = min(abs(s - s_crit));

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
        w0 = [k0;zeros(2,N)];

        if length(w_init) ~= 3
            w_init = w0(:,1);
        end

        k0prime = zeros(1,N);
        k0prime(1:ix_crit) = -2*kc1/L^2*(s_crit - s1) - 2*kc2/L*(1 - s_crit/L);
        k0prime(ix_crit+1:N) = -2*kc2/L*(1 - s2/L);
        w0prime = [k0prime;zeros(2,N)];

        %% integrate to calculate omega vector
        if nargout == 1
            wv = fn_intgEP_w0_Dimitri(w_init, w0, w0prime, 0, s0, ds, N, B, Binv);
        else
            [wv, pmat, Rmat] = fn_intgEP_w0_Dimitri(w_init, w0, w0prime, 0, s0, ds, N, B, Binv);
        end
    end
end