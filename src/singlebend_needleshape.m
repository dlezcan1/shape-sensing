%% singlebend_needleshape
% this is a function to produce the singlebend needle shape
% given the measured curvatures of the needle
%
% - written by: Dimitri Lezcano

function pos = singlebend_needleshape(curvatures, aa_tip_locs, L, kc_i, w_init_i, theta0)
% Input:
%   - curvatures: list of x-y curvatures measured at each of the AA locations
%           ( a #AA x 2 matrix ) (1/m)
%   - aa_locs: list of the active area locations (measured from the tip of
%   the needle) corresponding w/ rows of curvatures
%   - L: the needle length (mm)
%   - kc_i: the initial guess of kappa_c (double)
%   - w_init_i: the initial guess of w_init (3 x 1 vector, default is [kc_init; 0; 0])
    
    %% Arguments
    arguments
        curvatures (:, 2) {mustBeNumeric};
        aa_tip_locs (1,:) {mustBeNumeric};
        L double;
        kc_i double;
        w_init_i (3, 1) {mustBeNumeric} = [kc_i; 0; 0]; 
        theta0 double = 0;
    end

    %% material properties
    % Stainless Steel 304
    Emod = 200e9*1e-6; % 200 GPa, conversion from N/m^2 to N/mm^2
    Pratio = 0.29; % Poisson's ratio
    diam = 0.9; % in mm
    Ibend = pi*diam^4/64;

    Gmod = Emod/2/(1+Pratio);
    Jtor = pi*diam^4/32;

    BendStiff = Emod*Ibend;
    TorStiff = Gmod*Jtor;

    B = diag([BendStiff,BendStiff,TorStiff]);
    Binv = inv(B);

    
    %% Needle arclength set-up
    aa_base_locs = L - aa_tip_locs;
    ds = 0.5;
    s = 0:ds:L;
    
    % get the arclength indices that are valid
    aa_base_locs_valid = aa_base_locs(aa_base_locs >= 0);
    [~, s_idx_aa] = min(abs(s' - aa_base_locs_valid));
    curvs_aa = curvatures(aa_base_locs >= 0, :)*1e-3; % convert curvatures to 1/mm
    s_aa = s(s_idx_aa);
        
    
    %% Determine w_init and kc from measured curvatures (optimization)
    % intrinsic curvature values
    k0 = kc_i * (1 - s/L).^2;
    w0 = [k0; zeros(2, length(s))];
    
    k0_prime = -2*kc_i/L * (1 - s/L);
    w0_prime = [k0_prime; zeros(2, length(s))];
    
    % initial cost values
    eta = [w_init_i; kc_i];
%     wv = fn_intgEP_v1(w_init_i, w0, w0_prime, 0, ds, length(s), B, Binv);
    scalef0 = 1;
    Cval = costfn_shape_singlebend(eta, curvs_aa', s_idx_aa, ds, length(s), B, Binv, scalef0);
    scalef = 1/Cval;
    
    % optimization
    x0 = eta; % initial value
    LB = [-0.01*ones(3,1);0]; % lower bound
    UB = [0.01*ones(3,1);0.01]; % upper bound
    
    oldopts = optimset('fmincon');
    options = optimset(oldopts,'Algorithm','interior-point','TolFun',1e-6,'TolX',1e-8,...
        'MaxFunEvals',10000);
    [x, fval, exitflag] = fmincon( @(x) costfn_shape_singlebend(x, curvs_aa',...
        s_idx_aa, ds, length(s), B, Binv, scalef),...
        x0, [], [], [], [], LB, UB, [], options);
    
    % unpack optimization results
    w_init = x(1:3);
    kc = x(4);
    
    %% Generate needle shape
    w0 = kc * (1 - s/L).^2;
    w0_prime = -2*kc/L * (1 - s/L);
    
    
    
    
end
    