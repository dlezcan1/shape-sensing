function cost = cost_fn_kc_analysis_layer_v2(B, Binv, kc_init1, kc_init2, w_init, L_init, lengths, z_crit, p1, p2)
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

idx_proj = find( lengths ~= L_init);
idx_ref = find(lengths == L_init);
assert(length(idx_ref) == 1, "Length of ref. idx should be 1.");

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
    cost = mean(vecnorm(deviation_metric_data(positions(idx_proj), positions{idx_ref})));
    
end



    
    
    