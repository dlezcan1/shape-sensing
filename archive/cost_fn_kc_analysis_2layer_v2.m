function cost = cost_fn_kc_analysis_2layer_v2(B, Binv, kc_init1, kc_init2, w_init, L_init, lengths, z_crit, p1, p2)
% using deviation metric
% function to generate the cost for the 3 length scenario
% for the kappa_c analysis
% the third length will be 90 mm, which will be read in
% 
% kc_init1 - The 90mm kappa_c1 layer 
% kc_init2 - the 90mm kappa_c2 layer 
% w_init1  - the initial rotation insertion
% L1       - the 90 mm insertion length
% lengths  - the lengths to be optimized
% z_crit   - the critical z-value for the 2nd layer
% p1       - the exponent of the length term in 1st layer
% p2       - the exponent of the length term in 2nd layer
% 
% - written by: Dimitri Lezcano

%% Instantiations
    ds = 0.5;
    s = 0:ds:max(lengths);
    
    positions = cell(size(lengths));
    
    %% Get the Positions matrices for each of the lengths
    for i = 1:length(lengths)
    % wv instantations for EP equation
        kc_init1_i = kappa_c_p(kc_init1, L_init, lengths(i), p1);
        
        s_crit = determine_s_crit_1layer(kc_init1_i, w_init, z_crit, lengths(i), 0, ds, B, Binv);

    % get the position and add it to the cell
        positions{i} = position_kcp_2layer(B, Binv, kc_init1, kc_init2,...
            w_init, L_init, lengths(i), s_crit, p1, p2);
    
    end
    
%% Error Calculations
    point_cost =vecnorm(error_positions(positions));
    cost = mean(point_cost);
    
    
end


function cost = area_between_curves(pi, pj)
    % set-up for interpolation
    zi = pi(3,:);
    zj = pj(3,:);
    
    z_min = min(max(zi), max(zj));
    dz = 0.5;
    
    z = 0:dz:z_min;
    
    xi = interp1(zi, pi(1,:), z, 'spline');
    xj = interp1(zj, pj(1,:), z, 'spline');
    
    yi = interp1(zi, pi(2,:), z, 'spline');
    yj = interp1(zj, pj(2,:), z, 'spline');
    
    % take the area between curves.
    cost = sum(vecnorm([xi-xj; yi-yj]));
    
end    
    
    
    