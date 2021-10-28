%% singlebend_doublelayer_needleshape
% this is a function to produce the singlebend needle shape in 2-layer
% tissue given the measured curvatures of the needle.
%
% - written by: Dimitri Lezcano

function [pos, wv, Rmat, kc1, kc2, w_init, s_crit] = singlebend_doublelayer_needleshape(curvatures, aa_tip_locs, ...
                                needle_mechparams, L, z_crit, kc1_i, kc2_i, w_init_i, theta0, weights)
% Input:
%   - curvatures: list of x-y curvatures measured at each of the AA locations
%           ( a #AA x 2 matrix ) (1/m)
%   - aa_locs: list of the active area locations (measured from the tip of
%   the needle) corresponding w/ rows of curvatures
%   - L: the needle length (mm)
%   - z_crit: the length of the first tissue boundary
%   - kcX_i: the initial guess of kappa_c for layer X (double)
%   - w_init_i: the initial guess of w_init (3 x 1 vector, default is [kc_init; 0; 0])
%   - theta0: the offset insertion angle (double, default is 0)
    
    %% Arguments
    arguments
        curvatures (2, :) {mustBeNumeric};
        aa_tip_locs (1,:) {mustBeNumeric};
        needle_mechparams struct;
        L double;
        z_crit double;
        kc1_i double;
        kc2_i double = kc1_i;
        w_init_i (3, 1) {mustBeNumeric} = [kc_i; 0; 0]; 
        theta0 double = 0;
        weights (1,:) {mustBeEqualSize(weights,aa_tip_locs, 2)} = ones(1, length(aa_tip_locs));
    end

    %% material properties
    B = needle_mechparams.B;
    Binv = needle_mechparams.Binv;
    
    %% Needle arclength set-up
    aa_base_locs = L - aa_tip_locs;
    ds = 0.5;
    s = 0:ds:L;
    
    % get the arclength indices that are valid
    aa_base_locs_valid = aa_base_locs(aa_base_locs >= 0);
    [~, s_idx_aa] = min(abs(s' - aa_base_locs_valid));
    curvs_aa = curvatures(:, aa_base_locs >= 0)*1e-3; % convert curvatures to 1/mm
    curvs_aa = [curvs_aa; zeros(1,size(curvs_aa, 2))];
    weights = weights(aa_base_locs >= 0);
    s_aa = s(s_idx_aa);
        
    
    %% Determine w_init and kc from measured curvatures (optimization)
    % initial cost values
    eta = [w_init_i; kc1_i; kc2_i];
    scalef0 = 1;
    Cval = costfn_shape_singlebend_2layer(eta, curvs_aa, s_idx_aa, z_crit,...
        theta0, ds, length(s), B, Binv, scalef0, weights);
    scalef = 1/Cval;
    
    % optimization
    x0 = eta; % initial value
    LB = [-0.01*ones(3,1);0;0]; % lower bound
    UB = [0.01*ones(3,1);0.01;0.01]; % upper bound
    
    oldopts = optimset('fmincon');
    options = optimset(oldopts,'Algorithm','interior-point','TolFun',1e-8,'TolX',1e-8,...
        'MaxFunEvals',10000, 'Display', 'off');
    [x, fval, exitflag] = fmincon( @(x) costfn_shape_singlebend_2layer(x, curvs_aa,...
        s_idx_aa, z_crit, theta0, ds, length(s), B, Binv, scalef, weights),...
        x0, [], [], [], [], LB, UB, [], options);
    
    % unpack optimization results
    w_init = x(1:3);
    kc1 = x(4);
    kc2 = x(5);
    
    %% Generate needle shape
    [wv, pos, Rmat, s_crit] = fn_intgEP_zcrit_2layers_Dimitri(w_init, kc1, kc2, z_crit,...
        theta0, 0, ds, length(s), B, Binv);
    if max(pos(3,:)) <= z_crit
        kc2 = -1; 
    end
    
end

%% Helper functions
% cost function for needle shape
function y = costfn_shape_singlebend_2layer(eta,data,s_index_meas,z_crit, theta0, ds,N,B,Binv,scalef,weights) 
    weights = weights/sum(weights, 'all');
    % unpack the variables
    w_init = eta(1:3); 
    kc1 = eta(4); 
    kc2 = eta(5);

    % get the measured curvatures
    wv = fn_intgEP_zcrit_2layers_Dimitri(w_init, kc1, kc2, z_crit, theta0, 0, ds, N, B, Binv);
    
    % exclude torsion 
    yv = wv(1:2,s_index_meas) - data(1:2,:); 
    y = norm(yv.*weights,'fro')^2*scalef; 

end