%% main_Diffeq_FP_needle_erranalysis.m
%
% This script is to perform an error analysis on the shape sensing
%
% - written by: Dimitri Lezcano

clear

%%
global kc L w dw s_l w_i w_j w_k ds sigma A0 G0 

save_bool = true;

sigma = 0.001; % gaussian noise uncertainty

directory = "../data/";
file_base = directory + "DiffEq_Results_sigma-0.0025_ds-1.0_ideal";


%% Load the data
load(file_base + ".mat");


%% Get the mean omega
mean_w_FD = mean_prob_w(prob, w_i, w_j, w_k);

if isempty(w_init) % default w_init
    w_init = [kc; 0; 0];
end


%% create an ellipse matrix for each of the s values
prob_w12 = squeeze(sum(prob, 3)); % sum over w3 axis
% prob_w12 = prob_w12./sum(prob_w12, [1,2]); % normalize the probability 

[w1_grid, w2_grid] = meshgrid(w_i, w_j);

ellipse_mat = zeros(2,2,size(prob_w12,3));
mean_mat = zeros(2, size(prob_w12,3));
for l = 1:size(prob_w12, 3)
    % get the gaussian from 1 STD
    [mu, sigma_mat, ~] = fit_gaussian2d(w1_grid, w2_grid, prob_w12(:,:,l));
    mean_mat(:,l) = mu;
    ellipse_mat(:,:,l) = sigma_mat;
    
end

%% Get the angular deformations for the different parametriz
% theta_vals = 0:pi/4:2*pi-0.001; % the theta values we will use to parametrize the ellipse
theta_vals = linspace(0, 2*pi, 20); % the theta values we will use to parametrize the ellipse
u = [cos(theta_vals); sin(theta_vals)]; % unit vectors

w_err_mat = zeros(length(theta_vals), 3, length(ellipse_mat)); % #u x dim(w) x #arclength
w_err_mat(:,3,:) = repmat(reshape(mean_w_FD(3,:), 1, 1, []), length(theta_vals), 1); % assign torsion as mean
for i = 1:length(u) % iterate over unit vectors
    for l = 1:length(ellipse_mat) % iterate over all of the ellipses
        u_i = u(:,i); % the unit vector to work with
        sigma_mat = ellipse_mat(:,:,l); % the ellipse matrix
        
        % calculate the w value for w1 & w2
%         w12_i = sigma_mat * u_i + mean_mat(:,l);
        w12_i = sigma_mat * u_i + mean_w_FD(1:2,l);
    
        % append the result
        w_err_mat(i, 1:2, l) = w12_i;
        
    end
end

%% Start determining the shape integrations
% mean shape
w0 = mean_w_FD;
w0_prime = [diff(w0, 1, 2), zeros(3,1)];
[~,mean_shape,~] = fn_intgEP_w0_Dimitri(w_init, w0, w0_prime, 0, 0, ds, size(mean_mat, 2), B, Binv);


% shape bounds
shape_mat = zeros(size(w_err_mat));
for i = 1:size(w_err_mat, 1)
    % grab w0 values
    w0_i = squeeze(w_err_mat(i,:,:));
    
    % approximate w0_prime
    w0_prime_i = [diff(w0_i, 1, 2), zeros(3,1)];
    
    
    [~,pmat,~] = fn_intgEP_w0_Dimitri(w_init, w0_i, w0_prime_i, 0, 0, ds, size(shape_mat, 3), B, Binv);
    
    shape_mat(i,:,:) = pmat;
    
end

%% Plotting
% uncertainty shapes
for l = 1:size(shape_mat, 3)
    pts_l = shape_mat(:,:,l); % grab all the points in an arclength
    
    % plot the coordinates
    plot3(pts_l(:,3), pts_l(:,2), pts_l(:,1), 'k-', 'LineWidth', 2); hold on;
    
end
plot3(mean_shape(3,:), mean_shape(2,:), mean_shape(1,:),'r-');

hold off;
axis equal; grid on;
xlabel('z'); ylabel('y'); zlabel('x');
view([60, 15])


%% Saving
mean_w_FD = mean_mat;
mean_shape_FD = mean_shape;
sigma_w_mat_FD = ellipse_mat;
save('../data/FP-FD_ds-1.0-Shapes.mat', 'mean_shape_FD', 'sigma_w_mat_FD', 'mean_w_FD');