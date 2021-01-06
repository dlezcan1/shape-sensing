function [mean_error,z_interp, num_of_data_pts] = error_positions(positions)
% function to measure the errors along the needle
% takes the average of the positions at each of the z-values
%
% - written by: Dimitri Lezcano

% measure the mean error for each of the z-values

%% Instantiations
dz = 0.5;

%% Preparation of the loop data points
% determine the number of columns
z_tot_max = max([positions{1:end-1}],[],[1 2]); % assumes z-axis is longest
z_interp = 0:dz:z_tot_max;
N_cols = length(z_interp);

% set up the error data containers
mean_error = zeros(3,N_cols); % errors along the 2nd longest length
num_of_data_pts = zeros(1,N_cols);
%% loop through to get the errors
for i = 1:length(positions)
    for j = i+1:length(positions)
        posi = positions{i};
        posj = positions{j};
        
        z_max = min(max(posi(3,:)), max(posj(3,:)));
        
        % interpolation
        zij = z_interp(z_interp < z_max);
                
        xi = interp1(posi(3,:), posi(1,:), zij);
        xj = interp1(posj(3,:), posj(1,:), zij);
        
        yi = interp1(posi(3,:), posi(2,:), zij);
        yj = interp1(posj(3,:), posj(2,:), zij);
        
        % error processing
        error = [xi - xj; yi - yj; zij - zij]; % distance
        % sum error while casting to higher dimension
        mean_error = mean_error + [error zeros(3,N_cols - length(zij))];
        
        % count update of data points
        count_update = [ones(1,length(zij)) zeros(1,N_cols - length(zij))];
        num_of_data_pts = num_of_data_pts + count_update;
        
    end
end
% disp(N_cols);
% disp(num_of_data_pts);
mean_error = mean_error ./ num_of_data_pts;
