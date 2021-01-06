function cost = cost_fn_kc_analysis_v2(B, Binv, kc_init, w_init, L_init, lengths, p)
% function to generate the cost for the 3 length scenario
% for the kappa_c analysis
% the third length will be 90 mm, which will be read in
% 
% kc_init1 - The 90mm kappa_c 
% w_init1  - the initial rotation insertion
% L1       - the 90 mm insertion length
% L2, L3   - the two other insertion lengths of interest
% p        - the exponent of the length term
% 
% - written by: Dimitri Lezcano

%% Instantiations
    ds = 0.5;
    s = 0:ds:max(lengths);
    
    positions = cell(size(lengths));
    
    %% Get the Positions for each of the lengths
    for i = 1:length(lengths)
    % wv instantations for EP equation
        si = s(s <= lengths(i));

        kc_init_i = kappa_c_p(kc_init, L_init, lengths(i), p);
        w_init = [kc_init_i; 0; 0];

        k0_i = kc_init_i*(1 - si/lengths(i)).^2;
        k0_prime_i = -2*kc_init_i/lengths(i)*(1 - si/lengths(i));
        w0_i = [k0_i; zeros(2,length(si))];
        w0_prime_i = [k0_prime_i; zeros(2,length(si))];


    % EP calculations
        wv_i = fn_intg_EP_Dimitri(w_init, w0_i, w0_prime_i, lengths(i), si, ds, B, Binv);

        positions{i} = wv2r(wv_i, lengths(i));
    
    end
    
%% Error Calculations
    cost = mean(vecnorm(error_positions(positions)));
    
end


function cost = area_between_curves(pi, pj)
    % set-up for interpolation
    zi = pi(3,:);
    zj = pj(3,:);
    
    z_min = min(max(zi), max(zj));
    
    z = 0:.5:z_min;
    
    xi = interp1(zi, pi(1,:), z, 'spline');
    xj = interp1(zj, pj(1,:), z, 'spline');
    
    yi = interp1(zi, pi(2,:), z, 'spline');
    yj = interp1(zj, pj(2,:), z, 'spline');
    
    % take the area between curves.
    cost = sum(vecnorm([xi-xj; yi-yj]));
    
end
    
    
    
    