%% mean_prob_shape
% this script is to determine the mean probabilistic shape
%  based on the probability density of the omegas
%
% - written by: Dimitri Lezcano

directory = "Dimitri/Data/";
% data_file = directory + "DiffEq_Results_sigma_0.0001.mat";
data_file = directory + "DiffEq_Results_sigma_0.0005_bayes_meas.mat";
save_data_file = replace(data_file, ".mat", "_shapes.mat");

load(data_file); % load the data file

%% set-up for both shapes
ds_prob = ds;
BendStiff = Emod*Ibend;
TorStiff = Gmod*Jtor;

B = diag([BendStiff,BendStiff,TorStiff]);
Binv = inv(B);

if isempty(w_init)
    w_init = [kc; 0; 0];
    
end

%% find the mean probabilistic deformation
w_mean = mean_prob_w(prob, w_i, w_j, w_k);
w_mean_prime = [diff(w_mean, 1, 2), zeros(3,1)]; % derivative of mean omega


%% determine the expected and mean shape of the needle
% expected shape
[wv_exp, pmat_exp, ~] = fn_intgEP_v1_1layer(w_init,kc,0,0,ds,N,B,Binv);

% probabilistic shape
[wv_prob, pmat_prob, ~] = fn_intgEP_w0_Dimitri(w_init, w_mean, w_mean_prime,0,0,ds,N,B,Binv);

%% Plot the shapes
f3d = figure(1);
set(gcf,'units', 'normalized', 'position', [1/8, 1/4, 3/8, 3/8]);
f2d = figure(2);
set(gcf,'units', 'normalized', 'position', [1/2, 1/4, 3/8, 3/8]);

figure(1); % 3-D plots
plot3(pmat_exp(3,:), pmat_exp(1,:), pmat_exp(2,:),'linewidth',2, ...
        'DisplayName', 'expected'); hold on;
plot3(pmat_prob(3,:), pmat_prob(1,:), pmat_prob(2,:),'linewidth',2, ...
    'DisplayName', 'probabilistic'); hold off;
xlabel('z [mm]'); ylabel('x [mm]'); zlabel('y [mm]')
title('3-D view of Probabilistic Shape from F-P eq.');
legend();
grid on; axis equal;

figure(2); % 2-D plots
subplot(2,1,1);
plot(pmat_exp(3,:), pmat_exp(1,:),'linewidth',2, ...
   'DisplayName', 'expected'); hold on;
plot(pmat_prob(3,:), pmat_prob(1,:),'linewidth',2, ...
   'DisplayName', 'probabilistic'); hold off;
legend();
grid on; axis equal;

subplot(2,1,2);
plot(pmat_exp(3,:), pmat_exp(2,:),'linewidth',2, ...
   'DisplayName', 'expectd'); hold on;
plot(pmat_prob(3,:), pmat_prob(2,:),'linewidth',2, ...
   'DisplayName', 'expectd'); hold off;
xlabel('z [mm]'); ylabel('y [mm]');

sgtitle('2-D views of Probabilistic Shape from F-P eq.');
title('y-z axis view');
grid on; axis equal;

%% saving
savefig(f3d, replace(data_file, '.mat', '_shape_3d.fig'));
fprintf("Saved figure: %s\n", replace(data_file, '.mat', '_shape_3d.fig'));
saveas(f3d, replace(data_file, '.mat', '_shape_3d.png'));
fprintf("Saved image: %s\n", replace(data_file, '.mat', '_shape_3d.png'));
disp(' ');

savefig(f2d, replace(data_file, '.mat', '_shape_2d.fig'));
fprintf("Saved figure: %s\n", replace(data_file, '.mat', '_shape_2d.fig'));
saveas(f2d, replace(data_file, '.mat', '_shape_2d.png'));
fprintf("Saved image: %s\n", replace(data_file, '.mat', '_shape_2d.png'));