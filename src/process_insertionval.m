%% process_insertionval.m
%
% this is a script to run through the data points and generate the FBG shape
% from measurements
% 
% - written by: Dimitri Lezcano

%% Set-up
% directories to iterate through
expmt_dir = "../../data/needle_3CH_4AA_v2/Insertion_Experiment_04-22-21/";
trial_dirs = dir(expmt_dir + "Insertion*/");
mask = strcmp({trial_dirs.name},".") | strcmp({trial_dirs.name}, "..") | strcmp({trial_dirs.name}, "0");
trial_dirs = trial_dirs(~mask); % remove "." and ".." directories
trial_dirs = trial_dirs([trial_dirs.isdir]); % make sure all are directories

% files to find
fbgdata_file = "FBGdata.xls";
fbgdata_ref_wl_file = expmt_dir + "0_curvature/90_deg/" + fbgdata_file;
ref_wl_per_trial = true; % try not to use where possible

% weighted FBG measurement options
use_weights = true;

% saving options
save_bool = true;
fbgout_basefile = "FBGdata";
if use_weights == true
    fbgout_basefile = fbgout_basefile + "_FBG-weights";
end    
fbgout_posfile = fbgout_basefile + "_3d-position.xls";
fbgout_paramfile = fbgout_basefile + "_3d-params.txt";


% directory separation
if ispc
    dir_sep = '\';
else
    dir_sep = '/';
end

% calibraiton matrices file
calib_dir = "../../data/needle_3CH_4AA_v2/";
calib_file = calib_dir + "needle_params-Jig_Calibration_03-20-21_weighted.json";

% Initial guesses for kc and w_init
kc_i = 0.002;
w_init_i = [kc_i; 0; 0]; % ideal insertion
theta0 = 0;

%% Load the reference wavelengths
if ~ref_wl_per_trial
    ref_wls_mat = readmatrix(fbgdata_ref_wl_file, 'sheet', 'avg');
    ref_wls = mean(ref_wls_mat, 1); % the reference wavelengths
end

%% Load the calibration matrices and AA locations (from base)
fbgneedle = jsondecode(fileread(calib_file));

% AA parsing
num_aas = fbgneedle.x_ActiveAreas;
aa_base_locs_tot = struct2array(fbgneedle.SensorLocations); % total insertion length
weights = struct2array(fbgneedle.weights); 
aa_tip_locs = fbgneedle.length - aa_base_locs_tot;
cal_mats_cell = struct2cell(fbgneedle.CalibrationMatrices);
cal_mat_tensor = cat(3, cal_mats_cell{:});

if all(weights == 0)
    weights = ones(size(weights));
end

%% Iterate through the files
lshift = 1/6;
f3d = figure(1);
set(f3d,'units','normalized','position', [lshift + 0, 0.5, 1/3, .42])
f2d = figure(2);
set(f2d,'units','normalized','position', [lshift + 1/3, 0.5, 1/3, .42] )
f3d_insert = figure(3);
set(f3d_insert, 'units', 'normalized', 'position', [0, 0, 1/3, 0.42]);
f2d_insert = figure(4);
set(f2d_insert, 'units', 'normalized', 'position', [2/3, 0, 1/3, 0.42]); 
dir_prev = "";
for i = 1:length(trial_dirs)
    tic; 
    % trial operations
    L = str2double(trial_dirs(i).name);
    re_ret = regexp(trial_dirs(i).folder, "Insertion([0-9]+)", 'tokens');
    hole_num = str2double(re_ret{1}{1});
    
    % trial directory
    d = strcat(trial_dirs(i).folder,dir_sep, trial_dirs(i).name, dir_sep);
    fbg_file = d + fbgdata_file;
    
    % load the fbg shift in
    if ~ref_wl_per_trial
        wls_mat = readmatrix(fbg_file, 'Sheet', 'FBG_wavelength');
        wls_mat = wls_mat(all(wls_mat > 0, 2), :); % filter out any rows w/ 0 as FBG signal
        wls_mean = mean(wls_mat, 1);
        wls_shift = wls_mean - ref_wls;
        wls_shift = reshape(wls_shift, [], 3)'; % reshape the array so AA's are across rows and Ch down columns
    else
        wls_shift = readmatrix(fbg_file, 'Sheet', 'wavelength_shift');
        wls_shift = reshape(wls_shift, [], 3)';
    end
    
    
    % apply temperature compensation
    wl_shift_Tcorr = temperature_compensate(wls_shift);
    
    % use calibration senssors
    curvatures = calibrate_fbgsensors(wl_shift_Tcorr, cal_mat_tensor);
        
    % get the shape
    [pos, wv, Rmat, kc, w_init] = singlebend_needleshape(curvatures, aa_tip_locs, L, kc_i, w_init_i, theta0, weights);
    t = toc;
    
    % set new predictions
    kc_i = kc; 
    if i == 1
        w_init_i = w_init;
    elseif i > 1 && strcmp(trial_dirs(i).folder, trial_dirs(i-1).folder) 
        w_init_i = w_init;
    else
        w_init_i = [kc_i; 0; 0];
    end
    
    % plotting
    %- 3D
    figure(1);
    plot3(pos(3,:), pos(1,:), pos(2,:), 'linewidth', 2);
    axis equal; grid on;
    xlabel('z [mm]', 'FontWeight', 'bold'); ylabel('x [mm]', 'FontWeight', 'bold'); 
    zlabel('y [mm]', 'FontWeight', 'bold');
    title(d);
    
    %- 2D
    figure(2);
    subplot(2,1,1);
    plot(pos(3,:), pos(2,:), 'LineWidth', 2);
    xlabel('z [mm]', 'FontWeight', 'bold'); ylabel('y [mm]', 'FontWeight', 'bold');
    axis equal; grid on;
    
    subplot(2,1,2);
    plot(pos(3,:), pos(1,:), 'LineWidth', 2);
    xlabel('z [mm]', 'FontWeight', 'bold'); ylabel('x [mm]', 'FontWeight', 'bold');
    axis equal; grid on;
    
    % total insertion plots
    % - 3D total
    figure(3);
    if ~strcmp(dir_prev, trial_dirs(i).folder) % new trial
        hold off;
    end
    plot3(pos(3,:), pos(2,:), pos(1,:), 'linewidth', 2, 'DisplayName', sprintf("%.1f mm", L)); hold on;
    xlabel('z [mm]', 'FontWeight', 'bold'); ylabel('x [mm]', 'FontWeight', 'bold'); 
    zlabel('y [mm]', 'FontWeight', 'bold');
    legend();  
    axis equal; grid on;
    title(sprintf("Insertion #%d | FBG Shape Determination", hole_num));
    
    % - 3D total
    figure(4);
    subplot(2,1,1);
    if ~strcmp(dir_prev, trial_dirs(i).folder) % new trial
        hold off;
    end
    plot(pos(3,:), pos(2,:), 'LineWidth', 2, 'DisplayName', sprintf("%.1f mm", L)); hold on;
    xlabel('z [mm]', 'FontWeight', 'bold'); ylabel('x [mm]', 'FontWeight', 'bold');
    axis equal; grid on;
    
    subplot(2,1,2);
    if ~strcmp(dir_prev, trial_dirs(i).folder) % new trial
        hold off;
    end
    plot(pos(3,:), pos(1,:), 'LineWidth', 2, 'DisplayName', sprintf("%.1f mm", L)); hold on;
    xlabel('z [mm]', 'FontWeight', 'bold'); ylabel('x [mm]', 'FontWeight', 'bold');
    axis equal; grid on;
    legend()
    sgtitle(sprintf("Insertion #%d | FBG Shape Determination", hole_num));
    
    % save the data
    if save_bool
       % write position file
       writematrix(pos, d + fbgout_posfile);
       fprintf("Wrote 3D Positions: '%s'\n", d + fbgout_posfile);
       
       % write shape sensing parameters
       T = table(kc, w_init', 'VariableNames', {'kc', 'w_init'});
       writetable(T, d + fbgout_paramfile);
       fprintf("Wrote 3D Position Params: '%s'\n", d + fbgout_paramfile);
       
       % save figures
       fileout_base = strcat(trial_dirs(i).folder, dir_sep, fbgout_basefile);
       saveas(f3d_insert, fileout_base + "_3d-all-insertions.png");
       fprintf("Saved figure #%d: '%s'\n", f3d_insert.Number, ...
           fileout_base + "_3d-all-insertions.png");
       
       saveas(f2d_insert, fileout_base + "_2d-all-insertions.png");
       fprintf("Saved figure #%d: '%s'\n", f2d_insert.Number, ...
           fileout_base + "_2d-all-insertions.png");
       
    end
    
    % update previous directory
    dir_prev = trial_dirs(i).folder;
    
    % output
    fprintf("Finished trial: '%s' in %.2f secs.\n", d, t);
    disp(" ");
end

%% Completion
close all;
disp("Program Terminated.");

