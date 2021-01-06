function mean_error = deviation_metric_data(positions, pos_ref)
%% Preparation of the loop data points
    % determine the number of columns
%     z_tot_max = max([positions{1:end-1}],[],[1 2]); % assumes z-axis is longest
    dz = 0.5;
    z_tot_max = max(pos_ref(3,:));
    z_interp = 0:dz:z_tot_max;
    N_cols = length(z_interp);

    % set up the error data containers
    mean_error = zeros(3,N_cols); % errors along the 2nd longest length
    num_of_data_pts = zeros(1,N_cols);

%% loop through to get the errors
    for i = 1:length(positions)
        posi = positions{i};
     
        z_max = min(max(posi(3,:)), max(pos_ref(3,:)));
        
        % interpolation
        zir = z_interp(z_interp < z_max);
                
        xi = interp1(posi(3,:), posi(1,:), zir);
        xr = interp1(pos_ref(3,:), pos_ref(1,:), zir);
        
        yi = interp1(posi(3,:), posi(2,:), zir);
        yr = interp1(pos_ref(3,:), pos_ref(2,:), zir);

        error = [xi - xr; yi - yr; zir - zir];

        mean_error = mean_error + [error zeros(3,N_cols - length(zir))];
        
        count_update = [ones(1,length(zir)) zeros(1,N_cols - length(zir))];
        num_of_data_pts = num_of_data_pts + count_update;

end    

mean_error = mean_error ./ num_of_data_pts;