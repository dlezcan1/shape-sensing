%% main_Diffeq_FP_needle_gauss_meas.m
%
% Solve the F-P equation for the prob. density of deformation using Gaussian
% assumption using measurements
%
% - written by Dimitri Lezcano

%% Variable set-up
clear; clc;

%% saving options 
save_bool = true;

directory = "../data/";
file_base = directory + "DiffEq_Results_sigma-%.4f_gauss-approx_ds-%.1f_bayes-ideal";

%% preamble
% physical parameters
Emod = 200e9*1e-6; % 200 GPa, conversion from N/m^2 to N/mm^2
Pratio = 0.29; % Poisson's ratio
diam = 0.9; % in mm
Ibend = pi*diam^4/64;

Gmod = Emod/2/(1+Pratio);
Jtor = pi*diam^4/32;

A0 = Emod*Ibend; 
G0 = Gmod*Jtor; 

B = diag([A0, A0, G0]);
Binv = inv(B);

% insertion parameters
L = 90;
kc = 0.0025508; %0.003;
sigma = 10*0.001; % gaussian noise uncertainty
w_init = []; % ideal case
% w_init = [ 0.0035703; 0.00072161; -0.0086653 ]; % data insertion

if isempty(w_init)
    file_base = file_base + "_ideal";
else
    file_base = file_base + "_data";
end

% arclength coordinate
ds = 0.5;

% system parameters
%- constants
S.A0 = A0;
S.G0 = G0;
S.alpha1 = (S.A0 - S.G0)/S.A0;
S.alpha2 = - S.alpha1;
S.alpha3 = -S.A0 / S.G0;
S.B = B;
S.Binv = Binv;

%- insertion parameters
S.L = L;
S.kc = kc;
S.sigma = sigma; % standard deviation
S.sigma2 = S.sigma^2;
if isempty(w_init)
    S.w_init = [kc; 0; 0];
else
    S.w_init = w_init;
end
S.ds = ds;
S.arclengths = 0:S.ds:S.L;

%- measurement parameters
S.cov_m = [0.0520, 0.0532, inf]; % diagonal covariance of measurement probability
S.meas_locations = [10, 40, 75];
[~, S.meas_idxs] = min(abs(S.arclengths' - S.meas_locations), [], 1); % indexes just in case
S.appx_meas_locations = S.arclengths(S.meas_idxs); % approximate locations (if not already there)

%% Break up arclengths for each measurement cycle
S.arclengths_per_measurement = {};
idx_0 = 1;
for i = 1:length(S.meas_locations)
   idxs = idx_0:S.meas_idxs(i);
   S.arclengths_per_measurement{i} = S.arclengths(idxs);
   fprintf("AA: %d->%d | Min: %.2f | Max %.2f\n", i-1,i, min(S.arclengths_per_measurement{i}), ...
            max(S.arclengths_per_measurement{i}));
   idx_0 = S.meas_idxs(i) + 1;
end
S.arclengths_per_measurement{i+1} = S.arclengths(idx_0:end);
fprintf("AA: %d->tip | Min: %.2f | Max %.2f\n\n", i, min(S.arclengths_per_measurement{i}), ...
            max(S.arclengths_per_measurement{i}));

%% Dynamical formulation
needle_f = @(s, x) gaussian_dynamics(s, x, S);
x0 = [S.w_init; 0e-8*ones(3,1)];

% iterate through each of the measurements
x_init = x0;
s = []; x = [];
for i = 1:length(S.arclengths_per_measurement)
   % iterate up to next AA location
   [s_i, x_i] = ode45(@(s,x) needle_f(s,x), S.arclengths_per_measurement{i}, x_init);
   
   % Bayesian fusion from ideal sensor measurement
   if ismember(s_i(end), S.appx_meas_locations)
        xp_AAi = posterior_measurement_ideal(s_i(end), x_i(end,:), S);
        x_i(end,:) = reshape(xp_AAi, 1, []);
   end
   
   % append new data
   s = [s; s_i];
   x = [x; x_i];
   
   % update x_init
   x_init = reshape(x_i(end,:), [], 1);
end

mu = x(:,1:3)';
sigma2_vect = x(:,4:6)';
Sigma = zeros(3,3,size(mu, 2));
for i = 1:size(Sigma, 3)
   Sigma(:,:,i) = diag(sigma2_vect(:, i)); 
end

%% Generate the probabilistic needle shapes
% mean shape
w0 = mu;
w0_prime = [diff(w0, 1, 2), zeros(3,1)];
[~,mean_shape,~] = fn_intgEP_w0_Dimitri(S.w_init, w0, w0_prime, 0, 0, ds, size(mu, 2), B, Binv);

% sample Gaussian omega (w3 is assumed to be fixed)
theta_vals = linspace(0, 2*pi-0.001, 50);
u = [cos(theta_vals); sin(theta_vals)];
w_err_mat = zeros(length(theta_vals), 3, size(Sigma, 3));

for i = 1:length(theta_vals)
    u_i = u(:,i);
    for l = 1:size(Sigma, 3)
        sig_i = Sigma(1:2,1:2,l); % only grab first two
        
        w12_i = sig_i * u_i + mu(1:2,l);
        
        w_err_mat(i, 1:3, l) = [w12_i; mu(3,l)];
        
    end
end 

% shape bounds
shape_mat = zeros(size(w_err_mat));
for i = 1:size(w_err_mat, 1)
    % grab w0 values
    w0_i = squeeze(w_err_mat(i,:,:));
    
    % approximate w0_prime
    w0_prime_i = [diff(w0_i, 1, 2), zeros(3,1)];
    
    % generate sampled shape
    [~,pmat,~] = fn_intgEP_w0_Dimitri(S.w_init, w0_i, w0_prime_i, 0, 0, ds, size(shape_mat, 3), B, Binv);
    
    shape_mat(i,:,:) = pmat;
    
end

%% save the current run
if save_bool
    file_base = sprintf(file_base, sigma, S.ds);
    save(file_base + ".mat");
    disp("Saved: " + file_base + ".mat")
end

%% Plotting
% uncertainty shapes
fig_ws = figure(1);
for l = 1:size(shape_mat, 3)
    pts_l = shape_mat(:,:,l); % grab all the points in an arclength
    
    % plot the coordinates
    plot3(pts_l(:,3), pts_l(:,2), pts_l(:,1), 'k-', 'LineWidth', 2); hold on;
        
end
plot3(mean_shape(3,:), mean_shape(2,:), mean_shape(1,:),'r-'); hold off;
axis equal; grid on;
xlabel('z'); ylabel('y'); zlabel('x');
view([60, 15])

% 2-D overlay
fig_2d = figure(2);
for l = 1:size(shape_mat,3)
    pts_l = shape_mat(:,:,l);
    
    subplot(2,1,1);
    plot(pts_l(:,3), pts_l(:,1), 'k.'); hold on;
    
    subplot(2,1,2);
    plot(pts_l(:,3), pts_l(:,2), 'k.'); hold on;
    
end

subplot(2,1,1);
plot(mean_shape(3,:), mean_shape(1,:), 'r-'); hold off;
ylabel('x [mm]');
grid on; 

subplot(2,1,2);
plot(mean_shape(3,:), mean_shape(2,:), 'r-'); hold off;
xlabel('z [mm]'); ylabel('y [mm]');
grid on; 
sgtitle('2D view of 3-D Probailistic shape');

% cov in 3-D positions
dshape_mat = shape_mat - reshape(mean_shape, 1, size(mean_shape,1), size(mean_shape,2));
dshape_mat_norm = vecnorm(dshape_mat(:,[1;3],:), 2, 2).*sign(dshape_mat(:,1,:)); % out-of-plane
dshape_mat_norm(:,2,:) = vecnorm(dshape_mat(:,[2;3],:), 2, 2).*sign(dshape_mat(:, 2, :)); % in-plane
max_dshape_mat_norm = squeeze(max(dshape_mat_norm, [], 1));
min_dshape_mat_norm = squeeze(min(dshape_mat_norm, [], 1));

fig_cov = figure(3);
ax1 = subplot(2,1,1);
%- out-of-plane
plot(s, max_dshape_mat_norm(1,:), 'g-', 'LineWidth', 2, 'DisplayName', 'out-of-plane'); hold on;
plot(s, min_dshape_mat_norm(1,:), 'g-', 'LineWidth', 2, 'DisplayName', ''); hold on;

%- in-plane
plot(s, max_dshape_mat_norm(2,:), 'b--', 'LineWidth', 2, 'DisplayName', 'in-plane'); hold on;
plot(s, min_dshape_mat_norm(2,:), 'b--', 'LineWidth', 2, 'DisplayName', ''); hold on;

for i = 1:length(S.meas_locations)
    xline(S.meas_locations(i), 'r--', 'Linewidth', 2, 'DisplayName', "AA" + i); hold on;
end
title("Positional Variance against Arclength");
ylabel('Variance Bounds [mm]');
grid on; legend('Location', 'nw'); hold off;

% covariance from w
ax2 = subplot(2,1,2);
plot(s, sigma2_vect(1,:), '-', 'LineWidth', 2, 'DisplayName', '\omega_1'); hold on;
plot(s, sigma2_vect(2,:), '--', 'LineWidth',2, 'DisplayName',  '\omega_2'); hold on;
plot(s, sigma2_vect(3,:), '-', 'LineWidth',2, 'DisplayName',  '\omega_3'); hold on;

for i = 1:length(S.meas_locations)
    xline(S.meas_locations(i), 'r--', 'Linewidth', 2, 'DisplayName', "AA" + i); hold on;
end
hold off;
title("\omega Covariance");legend('Location', 'nw'); 
xlabel('s [mm]'); ylabel('\sigma^2_\omega');
grid on;

%% Saving
mean_w_gauss = mu;
mean_shape_gauss = mean_shape;
sigma_w_mat_gauss = Sigma;
save(sprintf('../data/FP-Gauss_ds-%.2f-Shapes.mat', S.ds), 'mean_shape_gauss', 'sigma_w_mat_gauss', 'mean_w_gauss');

saveas(fig_ws, file_base + '_ws-shape-3d.png');
disp("Saved figure: " +  file_base + '_ws-shape-3d.png');

saveas(fig_2d, file_base + '_ws-shape-2d.png');
disp("Saved figure: " +  file_base + '_ws-shape-2d.png');

saveas(fig_cov, file_base + '_cov.png');
disp("Saved figure: " +  file_base + '_cov.png');

%% Functions
% kappa_0 functions
function k0 = kappa_0(s, S)
    k0 = S.kc * (1 - s/S.L).^2;
    
end

function k0_prime = kappa_0_prime(s, S)
    k0_prime = -2 * S.kc/S.L * (1 - s/S.L);

end

% mean and covariance dynamics
function dx = gaussian_dynamics(s, x, S)
% Args:
%   s = the arclength currently at
%   x = (mu; sigma^2) where mu and sigma^2 are the current dynamics

    % unpack x
    mu = x(1:3);
    sigma2 = x(4:6);
    
    % prepare dynamics
    k0 = kappa_0(s, S);
    k0_prime = kappa_0_prime(s, S);
    
    % dynamics formulas
    dmu = [k0_prime + S.alpha1 * mu(2) * mu(3);
           k0 * mu(3) + S.alpha2 * mu(1) * mu(3);
           S.alpha3 * k0 * mu(2)];
       
    dsigma2 = S.sigma2 * ones(3,1);
    
    % concatenate dynamics
    dx = [dmu; dsigma2];

end

% calculate posterior measurement update based on ideal deformation
function [xp, mu_p, cov_p] = posterior_measurement_ideal(s, x, S)
    % unpack x
    mu = reshape(x(1:3), 3,1);
    cov = reshape(x(4:6), 3,1);

    % generate measurement mean
    mu_m = kappa_0(s, S) * [1; 0; 0];
    cov_m = S.cov_m;

    [mu_p, cov_p] = gauss_bayes_diagonal(mu, cov, mu_m, cov_m);
    xp = [mu_p; cov_p];

end