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
    
    cost = cost_shape_error(mean_shape, bound_shapes, 'type', params.type);

end

%% Helper functions
% shape error measurment
function cost = cost_shape_error(mean_shape, bound_shapes, kwargs)
    arguments
        mean_shape (3,:);
        bound_shapes (:,3,:) {mustBeEqualSize(mean_shape, bound_shapes, [2,3])};
        kwargs.type string = "tip-mean";
    end
    
    % Tip error 
    if strcmp(kwargs.type, "tip-mean")
        cost = tip_error_mean(mean_shape, bound_shapes);
    
    elseif strcmp(kwargs.type, "tip-bounds")
        cost = tip_error_bounds(mean_shape, bound_shapes);
        
    else
        error("NotImplementedError: '%s' is not a supported error type", kwargs.type);
    
    end
end

% center offset error
function err = tip_error_mean(mean_shape, bound_shapes)
    tip_bounds = bound_shapes(:,:,end)'; % N x 3 -> 3 x N
    tip_mean = mean_shape(:,end); % 3 x 1;
    
    % fit an ellipse to the bound points
    c = fit_ellipse2d(tip_bounds);
    
    err_vect = tip_mean - c;
    
    err = dot(err_vect, err_vect);
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
    integrand = @(t) dot(tip_mean - Sig*u(t) - c, tip_mean - Sig*u(t) - c);
    
    err = integral(@(t) integrand(t), 0, 2*pi);
    
end