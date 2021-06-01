%% FBGNeedle.m
% a class for FBGNeedle parametrizations
%
% not needed, use the py class to import!
%
% - written by: Dimitri Lezcano

classdef FBGNeedle
    properties
        length double               % length of needle
        num_channels int32          % number of channels
        num_aa int32                % number of AA's
        sensor_locations struct     % the sensor locations measured from the base
        calibration_matrices struct % associated calibration mats
        
    end
    
    
end