%% doublebend_singlelayer_needleshape
% this is a function to produce the singlebend needle shape
% given the measured curvatures of the needle
%
% - written by: Dimitri Lezcano

function [pos, wv, Rmat, kc, w_init] = doublebend_singlelayer_needleshape(curvatures, aa_tip_locs, ...
                                needle_mechparams, s_turn, L, kc_i, w_init_i, theta0, weights)
% Input:
%   - curvatures: list of x-y curvatures measured at each of the AA locations
%           ( a #AA x 2 matrix ) (1/m)
%   - aa_locs: list of the active area locations (measured from the tip of
%   the needle) corresponding w/ rows of curvatures
%   - s_turn: the length that the needle was turned.
%   - L: the needle length (mm)
%   - kc_i: the initial guess of kappa_c (double)
%   - w_init_i: the initial guess of w_init (3 x 1 vector, default is [kc_init; 0; 0])
%   - theta0: the offset insertion angle (double, default is 0)
    
    %% Arguments
    arguments
        curvatures (2, :) {mustBeNumeric};
        aa_tip_locs (1,:) {mustBeNumeric};
        needle_mechparams struct;
        s_turn double;
        L double;
        kc_i double;
        w_init_i (3, 1) {mustBeNumeric} = [kc_i; 0; 0]; 
        theta0 double = 0;
        weights (1,:) {mustBeEqualSize(weights,aa_tip_locs, 2)} = ones(1, length(aa_tip_locs));
    end
    
    %% Argument validation
    if s_turn <= L
        error("turn value is <= L, use single-bend single-layer function.");
    end

    %% material properties
    B = needle_mechparams.B;
    Binv = needle_mechparams.Binv;
    
    %% Needle arclength set-up
    aa_base_locs = L - aa_tip_locs;
    ds = 0.5;
    s = 0:ds:L;
    
    % determine the turn location
    [s_turn_rnd, s_idx_turn] = min(abs(s - s_turn));
    
    % get the arclength indices that are valid
    aa_base_locs_valid = aa_base_locs(aa_base_locs >= 0);
    [~, s_idx_aa] = min(abs(s' - aa_base_locs_valid));
    curvs_aa = curvatures(:, aa_base_locs >= 0)*1e-3; % convert curvatures to 1/mm
    curvs_aa = [curvs_aa; zeros(1,size(curvs_aa, 2))];
    weights = weights(aa_base_locs >= 0);
    s_aa = s(s_idx_aa);
        
    
    %% Determine w_init and kc from measured curvatures (optimization)
    % initial cost values
    eta = [w_init_i; kc_i];
    scalef0 = 1;
    Cval = costfn_shape_doublebend_1layer(eta, curvs_aa, s_idx_aa, s_idx_turn, ds, length(s), B, Binv, scalef0, weights);
    scalef = 1/Cval;
    
    % optimization
    x0 = eta; % initial value
    LB = [-0.01*ones(3,1);0]; % lower bound
    UB = [0.01*ones(3,1);0.01]; % upper bound
    
    oldopts = optimset('fmincon');
    options = optimset(oldopts,'Algorithm','interior-point','TolFun',1e-8,'TolX',1e-8,...
        'MaxFunEvals',10000, 'Display', 'off');
    [x, fval, exitflag] = fmincon( @(x) costfn_shape_doublebend_1layer(x, curvs_aa,...
        s_idx_aa, ds, length(s), B, Binv, scalef, weights),...
        x0, [], [], [], [], LB, UB, [], options);
    
    % unpack optimization results
    w_init = x(1:3);
    kc = x(4);
    
    %% Generate needle shape
    w0 = kc * (1 - s/L).^2;
    w0_prime = -2*kc/L * (1 - s/L);
    
    [wv, pos, Rmat] = fn_intgEP_w0_Dimitri(w_init, w0, w0_prime, theta0, 0, ds, length(s), B, Binv);
    
    
end

%% Helper functions
% cost function for needle shape
function y = costfn_shape_doublebend_1layer(eta,data,s_index_meas,s_idx_turn,ds,N,B,Binv,scalef, weights)
    weights = weights/sum(weights, 'all');
        
    w_init = eta(1:3);
    kc = eta(4);

    L = (N-1)*ds; % in mm
    s = 0:ds:L;

    s1 = s(1:s_idx_turn);
    s2 = s(s_idx_turn:end);

    kc1 = kc*((s1(end) - s1(1))/L)^(2/3);
    kc2 = kc*((s2(end) - s2(1))/L)^(2/3);

    % intrinsic curvature kappa_0 (quadratic)
    k0_1 = kc1*(1 - s1/L).^2;
    k0_2 = -kc2*(1 - s2/L).^2;
    k0_turn = (k0_1(end) + k0_2(1))/2; %0;
    k0 = [k0_1(1:end-1),k0_turn,k0_2(2:end)];

    k0prime1 = -2*kc1/L*(1 - s1/L);
    k0prime2 = 2*kc2/L*(1 - s2/L);
    k0prime_peak = (k0_2(2) - k0_1(end-1))/2/ds; 
    k0prime = [k0prime1(1:end-1),k0prime_peak,k0prime2(2:end)];

    % intrinsic curvature \omega_0
    w0 = [k0;zeros(size(s));zeros(size(s))];
    w0prime = [k0prime;zeros(size(s));zeros(size(s))];
    wv = fn_intgEP_w0_Dimitri(w_init,w0,w0prime,ds,N,B,Binv);

    % exclude torsion
    yv = wv(1:2,s_index_meas) - data(1:2,:);
    y = norm(yv.*weights,'fro')^2*scalef;


end