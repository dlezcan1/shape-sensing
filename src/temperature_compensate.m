%% tempatature_compensate.m
%
% this is a function to perform temperature compensation on FBG sensor readings
%
% - written by: Dimitri Lezcano

%% Main Function
function fbg_shift_Tcorr = temperature_compensate(fbg_shift)
% The # channels is assumed to be 3 at the moment 
%
% Input:
%   - fbg_signals: the fbg signals ( # channels x # AAs )
%   - cal_matrices: the calibration matrices per AA ( # channels x 2 x # AAs )
%
% Return:
%   - curvatures: the calibrated curvatures for [ 2 x #AAs ]
    %% Arguments
    arguments
        fbg_shift (3,:) {mustBeNumeric};
    end
    
    %% Temperature compensation
    fbg_shift_Tcorr = fbg_shift - mean(fbg_shift, 2); % subtract out mean of AA's
    
end