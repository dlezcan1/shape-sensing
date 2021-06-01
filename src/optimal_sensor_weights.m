%% optimal_sensor_weights.m
%
% this is a script to determine the optimal sensor weighting for the 
% FBG sensors by solving a least squares problem
%
% - written by: Dimitri Lezcano

%% Set-up 
% options
save_bool = true;

% python set-up
if ispc % windows file system
    pydir = "..\";
    
else
    pydir = "../";
    
end

if count(py.sys.path, pydir) == 0
    insert(py.sys.path, int32(0), pydir);
end

% file set-up
directory = "../../data/needle_3CH_4AA_v2/";
fbgneedle_param = directory + "needle_params-Jig_Calibration_03-20-21_weighted.json"; % weighted calibration
fbgneedle_param_weight = strrep(fbgneedle_param, '.json', '_weights.json'); % weighted fbg parmeters

% datadir = directory + "Jig_Calibration_08-05-20/"; % calibration data
datadir = directory + "Validation_Jig_Calibration_03-20-21/"; % validation data
data_mats_file = datadir + "Data Matrices_calval.xlsx";

if contains(fbgneedle_param, 'weighted')
    data_mats_proc_file = strrep(data_mats_file, '.xlsx', '_weighted.xlsx');
    fig_save_file = datadir + "CalWeight_Jig_Shape_fit";
else
    data_mats_proc_file = data_mats_file;
    fig_save_file = datadir + "Jig_Shape_fit";
end

data_mats_proc_file = strrep(data_mats_proc_file, '.xlsx', '_proc.xlsx');

% paramteter set-up
jig_offset = 26.0; % the jig offset of full insertion
% AA_weights = [];% [1, 0.9, 0.3, 0.0]; % [AA1, AA2, AA3, AA4] reliability weighting
% if ~isempty(AA_weights)
%     fig_save_file = fig_save_file + "_weighted";
%     
% end

%% Load FBGNeedle python class
fbg_needle = py.FBGNeedle.FBGNeedle.load_json(fbgneedle_param);
disp("FBGNeedle class loaded.")

% channel list
ch_list = 1:double(fbg_needle.num_channels);
CH_list = "CH" + ch_list;
aa_list = 1:double(fbg_needle.num_aa);
AA_list = "AA" + aa_list; % the "AAX" string version

%% load the processed data matrices
data_mats = struct();
for AA_i = AA_list
    data_mats.(AA_i) = readtable(data_mats_proc_file, 'Sheet', AA_i, ...
        'PreserveVariableNames', false, 'ReadRowNames', true); % remove the first column (exp #)
    disp(AA_i + " loaded.");
end
disp(' ');

%% Consolidate the readings into one W data matrix
% get the actual curvature we are trying to fit
W_act_aai = data_mats.(AA_i)(:, ["CurvatureX", "CurvatureY"]).Variables;
w_act = reshape(W_act_aai', [], 1);

% get the W data matrix
W = zeros(length(w_act), length(AA_list)); % the measurement array
for i = 1:length(AA_list)
    AA_i = AA_list(i);
    
    % get the data matrices for this AA
    W_meas_aai = data_mats.(AA_i)(:, ["PredCurvX", "PredCurvY"]).Variables;
    
    % add it to the measurement array
    W(:, i) = reshape(W_meas_aai', [], 1);
    
    
end

% % append the 'sum(weights) = 1' constraint ( eta := weights)
% W = [W; ones(1, size(W, 2))];
% w_act = [w_act; 1];

%% perform the non-negative least squares fit
curv_weight = curv_weight_rule(abs(w_act),false);
D = sqrt(diag(curv_weight));
weights = lsqlin(D*W, D*w_act, [], [], ones(1, size(W, 2)), [1], [0, 0, 0.075, 0.025], []); % set lower bound
weights = weights./sum(weights); % normalize just in case

disp("Weights:");
disp([AA_list; reshape(weights, 1, [])] )
base_print = ['[ ', repmat('%f, ', 1, size(W, 2) - 1), '%f ]\n\n'];
fprintf(base_print, weights);

disp('Total error');
disp(norm((W * weights - w_act)))

%% save the weights to the fbg_needle 
py_dict_weights = py.dict();
for i=1:length(AA_list)
    py_dict_weights{AA_list(i)} = weights(i);
    
end
if save_bool
    fbg_needle.set_weights(py_dict_weights);
    fbg_needle.save_json(fbgneedle_param_weight);
    fprintf('Saved json parmater file: %s\n', fbgneedle_param_weight);
end
%% Functions
% function for performing weighted least squares in curvature
function curv_weight = curv_weight_rule(k, trivial)
    disp(nargin)
    curv_weight = ones(length(k), 1);
    
    if ~trivial 
        for i = 1:length(k)
            if abs(k(i)) <= 1
                curv_weight(i) = 1;

            else
                curv_weight(i) = 0.05;

            end
            
            if mod(i, 2) == 0
                angle = atan2(k(i), k(i-1));
                if rad2deg(angle) == 90 
                    curv_weight(i) = 1;
                    curv_weight(i-1) = 1;
                end
            end
        end
    end
end
