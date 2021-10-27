function [theta0, varargout] = solve_theta0(carm_shape, L, kc, w_init, lb, ub, ds)
    arguments
        carm_shape (3,:);
        L {mustBePositive};
        kc double;
        w_init (3,1) = [kc; 0; 0];
        lb = deg2rad(-15);
        ub = deg2rad(15);
        ds = 0.5;
    end
    %% material properties
    % Stainless Steel 304
    Emod = 200e9*1e-6; % 200 GPa, conversion from N/m^2 to N/mm^2
    Pratio = 0.29; % Poisson's ratio
    diam = 1.27; % in mm
    Ibend = pi*diam^4/64;

    Gmod = Emod/2/(1+Pratio);
    Jtor = pi*diam^4/32;

    BendStiff = Emod*Ibend;
    TorStiff = Gmod*Jtor;

    B = diag([BendStiff,BendStiff,TorStiff]);
    Binv = inv(B); 
    
    %% Setup optimization
    theta0_init = 0;
    scalef0 = 1;
    scalef = 1/costfn_shape_theta0(theta0_init, carm_shape, kc, w_init, B, Binv, scalef0, ds);
    
    %% optiization
    oldopts = optimset('fmincon');
    options = optimset(oldopts,'Algorithm','interior-point','TolFun',1e-8,'TolX',1e-8,...
        'MaxFunEvals',10000, 'Display', 'off');
    [theta0, fval, exitflag] = fmincon(@(theta) costfn_shape_theta0(theta, carm_shape,...
                                                    kc, w_init, B, Binv, scalef, ds, L),...
                                       theta0_init, [], [], [], [], lb, ub, [], options);
                                   
   %% Generate the optimized FBG needle shape
    s = 0:ds:L;
    w0 = kc * (1 - s/L).^2 .* [1;0;0];
    w0prime = -2*kc/L * (1 - s/L) .* [1;0;0];
    [wv, pos, Rmat] = fn_intgEP_w0_Dimitri(w_init, w0, w0prime, theta0, 0, ds, ...
                        length(s), B, Binv);
                    
    varargout{1} = wv; varargout{2} = pos; varargout{3} = Rmat;
    
    
    
end

%% Helper functions
function cost = costfn_shape_theta0(theta0, carm_shape, kc, w_init, B, Binv, scalef, ds, L)
    
    
    %% Generate FBG needle shape
    s = 0:ds:L;
    w0 = kc * (1 - s/L).^2 .* [1;0;0];
    w0prime = -2*kc/L * (1 - s/L) .* [1;0;0];
    [~, pos, ~] = fn_intgEP_w0_Dimitri(w_init, w0, w0prime, theta0, 0, ds, ...
                        length(s), B, Binv);
                    
    %% Compute errors
    cost_points = errorLength3D(carm_shape, pos);
    cost = scalef * mean(cost_points);
            
end

