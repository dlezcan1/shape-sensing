%% calibrate_fbgsensors.m
%
% this is a function for calibrating the sensors
%
% - written by: Dimitri Lezcano

function curvatures = calibrate_fbgsensors(fbg_signals, cal_matrices)
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
        fbg_signals (3,:) {mustBeNumeric};
        cal_matrices (2, 3, :) {mustBeNumeric, mustBeEqualSize(fbg_signals, cal_matrices, [2, 3]),...
                    mustBeEqualSize(fbg_signals, cal_matrices, [1, 2])};
    end
        
   %% Begin the calibration
   num_aa = size(fbg_signals, 2);
   curvatures = zeros(2, num_aa);
   
   for i = 1:num_aa
       signal = fbg_signals(:, i);
       cal_mat = cal_matrices(:,:,i);
       curvatures(:,i) = cal_mat * signal;
   end
    
end