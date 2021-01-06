function mean_error = error_s_positions(positions)
% function to measure the errors along the needle
% takes the average of the positions at each of the arclength
%
% This function assumes that the indexing starts and increments
% at the same position and size.
%
% This function also assumes that the positions cell array is 
% ordered in increasing insertion length order 
% (i.e. 60, 90, 120, etc.)
%
% - written by: Dimitri Lezcano

N_cols = size(positions{end-1},2);
mean_error = zeros(3, N_cols);
num_of_data_pts = zeros(3, N_cols);

for i = 1:length(positions)-1
    for j = 1+i:length(positions)
        posi = positions{i};
        posj = positions{j};
        
        N_max = min(size(posi,2), size(posj,2)); % max number of elements
        
        % error positions
        error = posi(:, 1:N_max) - posj(:, 1:N_max);
        
        % update sum error
        mean_error = mean_error + [error zeros(3, N_cols - N_max)];
        
        % count_update
        num_of_data_pts = num_of_data_pts + [ones(1, N_max) zeros(1, N_cols - N_max)];
        
    end
end

mean_error = mean_error ./ num_of_data_pts;

end

