%% Script for running the kappa_c vs L analysis
% - written by Dimitri Lezcano

clear all; close all;

directory = "Dimitri/Data/";

save = true;

%% material properties
% Stainless Steel 304
Emod = 200e9*1e-6; % 200 GPa, conversion from N/m^2 to N/mm^2
Pratio = 0.29; % Poisson's ratio
diam = 0.9; % in mm
Ibend = pi*diam^4/64;

Gmod = Emod/2/(1+Pratio);
Jtor = pi*diam^4/32;

BendStiff = Emod*Ibend;
TorStiff = Gmod*Jtor;

B = diag([BendStiff,BendStiff,TorStiff]);
Binv = inv(B);

%% Initializations
% w_init = [0.0048; 0; 0]; % taken from main singlebend run
N = 150;
dL = 5;
% kc_init = 0.00330; % taken from main singlebend run
kc_init = 0.0020;
% kc_init = 0.001;
% w_init = zeros(3,1); % initial insertion rotation
% w_init = [kc_init; 0; 0]; 
w_init = 0;

L_init = 90; % mm
lengths = [60, 90, 120, 150, 180]; % mm
% lengths = 60:dL:180; % mm

%% Model analysis 
% for p
low_bound = -1;
up_bound = 2;

Tol = 1e-14;
opts = optimset('fmincon');
options = optimset(opts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
        'MaxFunEvals',10000,'Display','notify');
    
cost_fn = @(p) cost_fn_kc_analysis_v2(B, Binv, kc_init, w_init, L_init, lengths, p);
p0 = 2/3;

[p, fval, exitflag, output] = fmincon(cost_fn, 1/3, [],[],[],[],low_bound,up_bound,[],options);

disp("p value, optimized")
disp(p)
% p = 0.65;
fprintf('Cost: %f\n', cost_fn(p));


%% Get the resulting positions
positions = cell(size(lengths));
for i = 1:length(positions)
    positions{i} = position_kcp(B, Binv, kc_init, w_init, L_init, lengths(i), p);
    
end

%% Plot the positions (3D)
f3d = figure();
for i = length(positions):-1:1
    pi = positions{i};
    plot3(pi(3,:),pi(1,:),pi(2,:),'Displayname',sprintf("L = %d mm", lengths(i))); hold on;

end

hold off;
title("Position Plots: 3D")
if length(lengths) < 7
    title("Position Plots: 2D")
    legend()
else
    title(sprintf("Position Plots (%d - %d mm | %d mm incs.): 3D", min(lengths), max(lengths), dL));
end

xlabel('z'); ylabel('x'); zlabel('z');
grid on;
axis equal;
view([-13.9832, 41.0820])

%% Plot the positions (2D)
f2d = figure();
for i = length(positions):-1:1
    pi = positions{i};
    plot(pi(3,:),pi(2,:),'Displayname',sprintf("L = %d mm", lengths(i))); hold on;
    
end

hold off;

if length(lengths) < 7
    title("Position Plots: 2D")
    legend()
else
    title(sprintf("Position Plots (%d - %d mm | %d mm incs.): 2D", min(lengths), max(lengths), dL));
end
xlabel('z'); ylabel('y');
grid on;
axis equal;

%% Plot the errors (2d)
ferr = figure('units','normalized','outerposition',[0 0 1 1]);
[mean_error, z_vals] = error_positions(positions);
norm_mean_error = vecnorm(mean_error);

fprintf("Mean of the Mean Errors: %.5f mm\n", mean(norm_mean_error));

plot(z_vals, norm_mean_error);

if length(lengths) < 7
    title(sprintf("Mean Error | Overall mean: %.4f",...
        mean(norm_mean_error)));
    
else
    title(sprintf("Mean Error (%d - %d mm | %d mm incs.) | Overall mean: %.4f",...
        min(lengths), max(lengths), dL, mean(norm_mean_error)));
end


xlabel('z'); ylabel('Error (mm)');
grid on;
% ylim([0, 1]);
ylim([0,1.1*max(norm_mean_error)]);

%% save the plots
if ~save
    return
end

filename_pval = sprintf("kc-optimization_%dlengths_kc_%.5f_p_%.6f",length(lengths),kc_init, p);

% save the figure plots
savefig(f2d, directory + filename_pval + "_2d.fig")
disp("Saved figure: "+ directory + filename_pval + "_2d.fig");

savefig(f3d, directory + filename_pval + "_3d.fig");
disp("Saved figure: "+ directory + filename_pval + "_3d.fig");

savefig(ferr, directory + filename_pval + "_err.fig")
disp("Saved figure: "+ directory + filename_pval + "_err.fig");

% save the figures as PNG images.
saveas(f2d, directory + filename_pval + '_2d.png');
disp("Saved image: "+ directory + filename_pval + "_2d.png");

saveas(f3d, directory + filename_pval + '_3d.png');
disp("Saved image: "+ directory + filename_pval + "_3d.png");

saveas(ferr, directory + filename_pval + "_err.png")
disp("Saved figure: "+ directory + filename_pval + "_err.png");