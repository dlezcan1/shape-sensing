%% cost_fn_winit_qoptim.m
%
% cost function for q optimization for w_init
%
% - written by: Dimitri Lezcano

function cost = cost_fn_winit_qoptim(q, kc_ref, winit_ref, L_ref, ...
                    needleparams, theta0, p, N_pred, dL, ds)
    %% Arguments Block
    arguments
        q double;
        kc_ref double;
        winit_ref (3,1);
        L_ref double;
        needleparams struct;
        theta0 double = 0;
        p double = 0.592;
        N_pred {mustBeInteger, mustBeGreaterThanOrEqual(N_pred, 5)} = 5;
        dL double {mustBePositive} = 5;
        ds double {mustBePositive} = 0.5;
    end
    
    %% Setup
    pmat_pred = cell(N_pred+1,1);
    lengths = L_ref + (0:dL:N_pred*dL);
    
    %% Generate needle shapes
    for i = 1:numel(lengths)
        [~, pmat_i, ~] = generate_predicted_insertion_1layer(lengths(i), L_ref, ...
            kc_ref, winit_ref, p, q, needleparams, 'ds', ds, 'theta0', theta0);
        pmat_pred{i} = pmat_i;
         
    end
    
    %% Generate cost of the function
    cost = TipAerror_Dimitri(pmat_pred, lengths);
    
end