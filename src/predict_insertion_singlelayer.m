%% predict_insertion_singlelayer.m
%
% this is a function to predict single layer insertion
%
% - written by: Dimitri Lezcano

function [p_pred, R_pred, wv_pred, q_best] = predict_insertion_singlelayer(L_pred, L_ref, ...
                                    kc_ref, w_init_ref, needleparams, p, params)
    arguments
        L_pred  {mustBePositive, mustBeNonempty};
        L_ref double;
        kc_ref double;
        w_init_ref (3,1);
        needleparams
        p double = 0.592;
        params.q0 double = p;
        params.dL double = 5
        params.ds double = 0.5;
        params.theta0 double = 0;
        params.optim_lb double = 0.0;
        params.optim_ub double = 1.0;
        params.optim_Tol double = 1e-14;
        params.optim_display string = 'notify';
    end
    
    %% Set-up
    % set up cells if there are more than one lengths
    if numel(L_pred) > 1
        p_pred = cell(size(L_pred));
        R_pred = cell(size(L_pred));
        wv_pred = cell(size(L_pred));
    end
    
    
    %% q param optimization for w_init
    N_pred = max(numel(L_pred), 5);
    
    % optimization cost function
    cost_fn = @(q) cost_fn_winit_qoptim(q, kc_ref, w_init_ref, L_ref, needleparams,...
                params.theta0, p, N_pred, params.dL, params.ds);
    
    % optimization options (fmincon)
    opts = optimset('fmincon');
    options = optimset(opts,'Algorithm','interior-point','TolFun',params.optim_Tol,...
            'TolX',1e-8,'MaxFunEvals',10000,'Display',params.optim_display);
        
    [q_best, fval, exitflag, output] = fmincon(cost_fn, params.q0, [], [], [], [], ...
            params.optim_lb, params.optim_ub, [], options);
    
    fmincon_results = struct('fval', fval, 'exitflag', exitflag, 'output', output);
    
    
    %% Prediction
    if numel(L_pred) == 1
        % generate predicted insertion
        [wv_pred, p_pred, R_pred] = generate_predicted_insertion_1layer(L_pred, L_ref, ...
                kc_ref, w_init_ref, p, q_best, needleparams, 'ds', params.ds, ...
                'theta0', params.theta0);
    else
        for i = 1:numel(L_pred)
            % generate predicted insertion
            [wv_pred_i, p_pred_i, R_pred_i] = generate_predicted_insertion_1layer(L_pred(i), L_ref, ...
                    kc_ref, w_init_ref, p, q_best, needleparams, 'ds', params.ds, ...
                    'theta0', params.theta0);
            
            % add new data
            wv_pred{i} = wv_pred_i;
            p_pred{i}  = p_pred_i;
            R_pred{i}  = R_pred_i;
        end
    end
end