%% parse_fbgneedle_json.m
%
% function to parse the fbgneedle json file for needle parameters
%
% - written by: Dimitri Lezcano

function fbgneedle = parse_fbgneedle_json(filename_json)
    arguments
        filename_json string;
    end
    
    % read in the Json file
    fbgneedle_json = jsondecode(fileread(filename_json));
    
    % Parse physical parameters
    fbgneedle.num_channels = fbgneedle_json.x_Channels;
    fbgneedle.num_activeAreas = fbgneedle_json.x_ActiveAreas;
    fbgneedle.length = fbgneedle_json.length;
    
    % parse AA sensor locations
    fbgneedle.slocs_base = struct2array(fbgneedle_json.SensorLocations);
    fbgneedle.slocs_tip = fbgneedle.length - fbgneedle.slocs_base;
    
    % parse AA reliability weightings
    fbgneedle.aa_weights = struct2array(fbgneedle_json.weights);
    if all(fbgneedle.aa_weights == 0)
        fbgneedle.aa_weights = ones(size(fbgneedle.aa_weights));
    end
    
    % parse calibration matrices
    cal_mats_cell = struct2cell(fbgneedle_json.CalibrationMatrices);
    fbgneedle.aa_calMats = cat(3, cal_mats_cell{:});
    
end
