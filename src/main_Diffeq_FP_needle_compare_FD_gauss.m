%% main_Diffeq_FP_needle_compare_FD_gauss.m
%
%  script to compare the results for the FD method to the Gaussian approximation
%  of the needle's Fokker-Planck equation
%
% - written by: Dimitri Lezcano

set(0,'DefaultAxesFontSize',24);

%% Set-Up
% file set-up
data_dir = "../data/";
fd_file = data_dir + "FP-FD-Shapes.mat";
gauss_file = data_dir + "FP-Gauss-Shapes.mat";
fileout_base = "compare_FD-Gauss";

% arclength parameters
ds = 1;
L = 90;
s = 0:ds:L;

%% Load the shapes
fd = load(fd_file);
gauss = load(gauss_file);


%% Compute errors
delta_shape = abs(fd.mean_shape_FD - gauss.mean_shape_gauss);
delta_sigma = abs(fd.sigma_mat_FD - gauss.sigma_mat_gauss);

% mean shape L2 norm
L2_shape = vecnorm(delta_shape);

% Frobenious norm for sigma
fro_sigma = zeros(1, length(delta_sigma));
for i = 1:length(fro_sigma)
    fro_sigma(i) = norm(delta_sigma(:,:,i), 'fro');
    
end

% L2 of eigenvalues
L2eigs_sigma = zeros(1, length(delta_sigma));
delta_eigs = zeros(2, length(delta_sigma));
for i = 1:length(L2eigs_sigma)
    eig_fd = eig(fd.sigma_mat_FD(:,:,i));
    eig_gauss = eig(gauss.sigma_mat_gauss(:,:,i));
    delta_eigs(:,i) = abs(eig_fd - eig_gauss);
    L2eigs_sigma(i) = norm(delta_eigs(:,i));
    
end

% determinant ratio for sigma (|Gauss|/|FD|)
detratio_sigma = zeros(1, length(delta_sigma));
for i = 1:length(detratio_sigma)
    detratio_sigma(i) = det(gauss.sigma_mat_gauss(:,:,i))/det(fd.sigma_mat_FD(:,:,i));
    
end

%% Generate probabilistic shapes
% Gaussian 
%% Generate the probabilistic needle shapes
% mean shape
w0 = gauss.mean_shape_gauss;
w0_prime = [diff(w0, 1, 2), zeros(3,1)];
[~,mean_shape,~] = fn_intgEP_w0_Dimitri(S.w_init, w0, w0_prime, 0, 0, ds, size(mu, 2), B, Binv);

% sample Gaussian omega (w3 is assumed to be fixed)
theta_vals = linspace(0, 2*pi-0.001, 50);
u = [cos(theta_vals); sin(theta_vals)];
w_err_mat = zeros(length(theta_vals), 3, size(Sigma, 3));

for i = 1:length(theta_vals)
    u_i = u(:,i);
    for l = 1:size(Sigma, 3)
        sig_i = 2*Sigma(1:2,1:2,l); % only grab first two
        
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

%% Display statistics
disp("Mean Shape Error Statistics");
fprintf(" Min: %f \n Mean: %f \n Max: %f \n\n", min(L2_shape), mean(L2_shape), max(L2_shape));

disp("Covariance Matrix Error Statistics");
disp(" Frobenius Norm");
fprintf("  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(fro_sigma), mean(fro_sigma), max(fro_sigma));
    
disp(" L2 Eigenvalues Norm");
fprintf("  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(L2eigs_sigma), mean(L2eigs_sigma), max(L2eigs_sigma));

disp(" |Cov_Gauss|/|Cov_FD|");
fprintf("  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(detratio_sigma), mean(detratio_sigma, 'double', 'omitnan'),...
        max(detratio_sigma));

%% Plotting
fig_err = figure(1);
subplot(1,2,1);
plot(s, delta_shape(1,:), 'DisplayName', 'x', 'LineWidth', 2); hold on;
plot(s, delta_shape(2,:), 'DisplayName', 'y', 'LineWidth', 2); hold on;
plot(s, delta_shape(3,:), 'DisplayName', 'z', 'LineWidth', 2); hold on;
% plot(s, L2_shape, 'DisplayName', 'L2', 'LineWidth', 2); hold on;
hold off;
grid on;
xlabel('arclength (mm)'); ylabel('Error (mm)');
legend('Location', 'best'); title('Mean Shape Errors');

subplot(1,2,2);
plot(s, fro_sigma, 'DisplayName', 'Fro', 'LineWidth', 2); hold on;
% plot(s, L2eigs_sigma, 'DisplayName', 'L2 of e-vals.', 'LineWidth', 2); hold on;
% plot(s, detratio_sigma, 'DisplayName', '|\Sigma_{G}|/|\Sigma_{FD}|', 'LineWidth', 2); hold on;
hold off;
grid on; 
xlabel('arclength (mm)'); ylabel('Error');
legend('Location', 'best'); title('Covariance Matrix Errors');

% 3-D workspace Gauss
fig_gauss = figure(2);

% 3-D workspace FD
fig_fd = figure(3);

 %% Saving
% Statistics 
fid = fopen(data_dir + fileout_base + "_stats.txt", 'w');
fprintf(fid, "Mean Shape Error Statistics\n");
fprintf(fid, " Min: %f \n Mean: %f \n Max: %f \n\n", min(L2_shape), mean(L2_shape), max(L2_shape));

fprintf(fid, "Covariance Matrix Error Statistics\n");
fprintf(fid, " Frobenius Norm\n");
fprintf(fid, "  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(fro_sigma), mean(fro_sigma), max(fro_sigma));
    
fprintf(fid, " L2 Eigenvalues Norm\n");
fprintf(fid, "  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(L2eigs_sigma), mean(L2eigs_sigma), max(L2eigs_sigma));

fprintf(fid, " |Cov_Gauss|/|Cov_FD|\n");
fprintf(fid, "  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(detratio_sigma), mean(detratio_sigma, 'double', 'omitnan'),...
        max(detratio_sigma));

fclose(fid);
fprintf("Saved statistics file: '%s'\n", data_dir + fileout_base + "_stats.txt");

% Error Figure
saveas(fig_err, data_dir + fileout_base + "_error.png");
fprintf("Saved Error figure: %s\n", data_dir + fileout_base + "_error.png");

savefig(fig_err, data_dir + fileout_base + "_error.fig");
fprintf("Saved Error figure: %s\n", data_dir + fileout_base + "_error.fig");

