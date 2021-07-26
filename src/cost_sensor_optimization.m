%% cost_sensor_optimization.m
%
% this is a function to provide a cost for the sensor optimization
%
% types available:
%   "tip-mean" - compares mean of tip and mean of gaussian
%   "tip-bounds" - compares error bounds from shape sensing
%
% - written by: Dimitri Lezcano

function cost = cost_sensor_optimization(slocs, params)
    %% Arguments block
    arguments
        slocs (1,:);
        params.kc double = 0.0025508;
        params.w_init = [];
        params.sigma double = 0.0025;
        params.L double = 90;
        params.ds double = 0.5;
        params.type string = "tip-mean";
    end
    
    %% Generate the needle shape
    [mean_shape, bound_shapes] = needle_gauss_meas(slocs, 'kc', params.kc, ...
        'w_init', params.w_init, 'sigma', params.sigma, 'L', params.L, ...
        'ds', params.ds);
    
    cost = cost_shape_error(mean_shape, bound_shapes, params);

end

%% Helper functions
% shape error measurment
function cost = cost_shape_error(mean_shape, bound_shapes, params)
    arguments
        mean_shape (3,:);
        bound_shapes (:,3,:) {mustBeEqualSize(mean_shape, bound_shapes, [2,3])};
        params struct;
    end
    
    % Tip error 
    
    switch params.type
        case "tip-mean"
            cost = tip_error_mean(mean_shape, bound_shapes);
            
        case "tip-bounds"
            cost = tip_error_bounds(mean_shape, bound_shapes);
        
        case "shape-bounds"
            cost = shape_error_bounds(mean_shape, bound_shapes, params);
        
        case "shape-mean"
            cost = shape_error_mean(mean_shape, bound_shapes, params);
            
        otherwise
            error("NotImplementedError: '%s' is not a supported error type", params.type);
    
    end
end

% volume error for all points
function err = shape_error_bounds(mean_shape, bound_shape, params)
    % preparations
    u = @(t) [cos(t); sin(t); zeros(1,numel(t))];
%     integrand = @(t, Sig, c, pt_mean) dot(pt_mean - Sig*u(t) - c, pt_mean - Sig*u(t) - c);
    integrand = @(t, Sig, c, pt_mean) vecnorm(pt_mean - Sig*u(t) - c);

    % run through all the points
    N = size(bound_shape, 3); % number of points
    err = 0; 
    for i = 1:N
        % find the fitting ellipse
        pt_mean_i = mean_shape(:,i);
        bounds_i = bound_shape(:,:,i)';
        [c, Sig2] = fit_ellipse2d(bounds_i);
        Sig2(isnan(Sig2)) = 0;
        Sig = sqrtm(Sig2);
        
        err = err + params.ds*integral(@(t) integrand(t, Sig, c, pt_mean_i), 0, 2*pi);
    end
end

function err = shape_error_mean(mean_shape, bound_shapes, params)
    N = size(mean_shape,2);
    err = 0;
    for i = 1:N
        pt_mean_i = mean_shape(:,i);
        bounds_i = bound_shapes(:,:,i)';
        c_i = fit_ellipse2d(bounds_i);
        
        err_vect = pt_mean_i - c_i;
        
        err = err + norm(err_vect) * params.ds;
    end

end

% area error
function err = tip_error_bounds(mean_shape, bound_shapes)
    tip_bounds = bound_shapes(:,:,end)';
    tip_mean = mean_shape(:,end);
    
    % fit an ellipse to the bound points
    [c, Sig2] = fit_ellipse2d(tip_bounds);
    Sig = sqrtm(Sig2);
    
    % Perform integration
    u = @(t) [cos(t); sin(t); zeros(size(t))];
%     integrand = @(t) dot(tip_mean - Sig*u(t) - c, tip_mean - Sig*u(t) - c);
    integrand = @(t) vecnorm(tip_mean - Sig*u(t) - c);
    
    err = integral(@(t) integrand(t), 0, 2*pi);
    
end

% center offset error
function err = tip_error_mean(mean_shape, bound_shapes)
    tip_bounds = bound_shapes(:,:,end)'; % N x 3 -> 3 x N
    tip_mean = mean_shape(:,end); % 3 x 1;
    
    % fit an ellipse to the bound points
    c = fit_ellipse2d(tip_bounds);
    
    err_vect = tip_mean - c;
    
    err = norm(err_vect);
end