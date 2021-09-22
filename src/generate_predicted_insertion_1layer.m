%% generate_predicted_insertion_1layer
% 
% this is a function to generate a predicted insertion in 1 layer tissue
% given reference parameters and p and q values
%
% - written by: Dimitri Lezcano

function [wv, pmat, Rmat] = generate_predicted_insertion_1layer(L, L_ref, kc_ref,...
        w_init_ref, p, q, needleparams, insertionparams)
    %% Arguments block
    arguments
        L double;
        L_ref double;
        kc_ref double;
        w_init_ref (3,1);
        p double;
        q double;
        needleparams struct;
        insertionparams.ds double {mustBePositive} = 0.5;
        insertionparams.theta0 double = 0;
    end
    
    %% Set-up
    % predicted parameters
    kc = kappa_c_p(kc_ref, L_ref, L, p);
    w_init = w_init_q(w_init_ref, L_ref, L, q);
    
    %% Generate the shape
    [wv, pmat, Rmat] = fn_intgEP_1layer_Dimitri(kc, w_init, insertionparams.theta0, ...
                            L, 0, insertionparams.ds, needleparams.B, needleparams.Binv);
                        
end