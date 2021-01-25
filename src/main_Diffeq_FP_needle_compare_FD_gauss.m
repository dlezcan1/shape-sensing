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
detratio_sigma(isnan(detratio_sigma)) = 0;

%% Display statistics
disp("Mean Shape Error Statistics");
fprintf(" Min: %f \n Mean: %f \n Max: %f \n\n", min(L2_shape), mean(L2_shape), max(L2_shape));

disp("Covariance Matrix Error Statistics");
disp(" Frobenius Norm");
fprintf("  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(fro_sigma), mean(fro_sigma), max(fro_sigma));
    
disp(" L2 Eigenvalues Norm");
fprintf("  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(L2eigs_sigma), mean(L2eigs_sigma), max(L2eigs_sigma));

disp(" |Cov_Gauss|/|Cov_FD|");
fprintf("  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(detratio_sigma), mean(detratio_sigma), max(detratio_sigma));

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
fprintf(fid, "  Min: %f \n  Mean: %f \n  Max: %f \n\n", min(detratio_sigma), mean(detratio_sigma), max(detratio_sigma));

fclose(fid);
fprintf("Saved statistics file: '%s'\n", data_dir + fileout_base + "_stats.txt");

% Error Figure
saveas(fig_err, data_dir + fileout_base + "_error.png");
fprintf("Saved Error figure: %s\n", data_dir + fileout_base + "_error.png");

savefig(fig_err, data_dir + fileout_base + "_error.fig");
fprintf("Saved Error figure: %s\n", data_dir + fileout_base + "_error.fig");

