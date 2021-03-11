%% main_Diffeq_FP_needle_FD_convergence.m
%
%  script to check for convergence of the results for the FD method 
%  of the needle's Fokker-Planck equation
%
% - written by: Dimitri Lezcano

set(0,'DefaultAxesFontSize',24);

%% Set-Up
% file set-up
data_dir = "../data/";
fd_ref_file = data_dir + "FP-FD_ds-0.5-Shapes.mat"; % refined ds = 0.5 mm
fd_file = data_dir + "FP-FD_ds-1.0-Shapes.mat"; % rough ds = 1.0
fileout_base = "FD_convergence";

% arclength parameters
ds = 1;
L = 90;
s = 0:ds:L;

%% Load the shapes
fd = load(fd_file);
fd_ref = load(fd_ref_file);


%% Compute errors
delta_shape = abs(fd.mean_shape_FD - fd_ref.mean_shape_FD(:,1:2:end));
delta_sigma = abs(fd.sigma_w_mat_FD - fd_ref.sigma_w_mat_FD(1:2,1:2,1:2:end));

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
    eig_fd = eig(fd.sigma_w_mat_FD(1:2,1:2,i));
    eig_fd_ref = eig(fd_ref.sigma_w_mat_FD(1:2,1:2,2*i-1));
    delta_eigs(:,i) = abs(eig_fd - eig_fd_ref);
    L2eigs_sigma(i) = norm(delta_eigs(:,i));
    
end

% determinant ratio for sigma (|FD|/|FD_ref|)
detratio_sigma = zeros(1, length(delta_sigma));
for i = 1:length(detratio_sigma)
    detratio_sigma(i) = det(fd.sigma_w_mat_FD(1:2,1:2,i))/det(fd_ref.sigma_w_mat_FD(:,:,2*i-1));
    
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

% - 2D mean shapes
f2d = figure(2);
subplot(2,1,1);
plot(fd.mean_shape_FD(3,:), fd.mean_shape_FD(1,:), '-', 'LineWidth', 2, 'DisplayName', 'ds = 1.0 mm'); hold on;
plot(fd_ref.mean_shape_FD(3,:), fd_ref.mean_shape_FD(1,:), '--','LineWidth', 2, 'DisplayName', 'ds = 0.5 mm'); hold off;
legend('Location', 'best');
ylabel('x [mm]');
axis equal; grid on;

subplot(2,1,2);
plot(fd.mean_shape_FD(3,:), fd.mean_shape_FD(2,:), '-', 'LineWidth', 2, 'DisplayName', 'ds = 1.0 mm'); hold on;
plot(fd_ref.mean_shape_FD(3,:), fd_ref.mean_shape_FD(2,:), '--','LineWidth', 2, 'DisplayName', 'ds = 0.5 mm'); hold off;
xlabel('z [mm]'); ylabel('y [mm]');
axis equal; grid on;
sgtitle('Mean Shapes for needle deformation: Convergence Test');


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

saveas(f2d, data_dir + fileout_base + "_2D-mean.png");
fprintf("Saved Mean Shape figure: %s\n", data_dir + fileout_base + "_2D-mean.png");

savefig(fig_err, data_dir + fileout_base + "_2D-mean.fig");
fprintf("Saved Mean Shape figure: %s\n", data_dir + fileout_base + "_2D-mean.fig");

