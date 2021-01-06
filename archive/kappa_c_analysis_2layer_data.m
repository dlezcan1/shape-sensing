%% Script for running the kappa_c vs L analysis
% - written by Dimitri Lezcano

clear all; close all;

directory = "Dimitri/Data/";
save = false;

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
N = 150;
dL = 1;
kc_init1 = 0.002;
kc_init2 = 0.003;
% kc_init2 = kc_init1;
z_crit = 45;
w_init = [ 0.0026392; -0.00016634;  0.0050071 ]; % initial boundary condition
p1 = 0.493;

L_init = 90; % mm
lengths = [60, 90, 120, 150]; % mm
% lengths = 60:dL:150; % mm

%% Model analysis 
% for p2
low_bound = -2;
up_bound = 2;

Tol = 1e-14;
opts = optimset('fmincon');
options = optimset(opts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
        'MaxFunEvals',10000,'Display','notify');
    

cost_fn = @(p) cost_fn_kc_analysis_2layer_data(B, Binv, kc_init1, kc_init2, w_init, L_init, lengths, z_crit, p1, p);
eta = p1;

[p2, fval, exitflag, output] = fmincon(cost_fn, eta, [],[],[],[],low_bound,up_bound,[],options);

disp("p2 value, optimized:")
disp(p2)
% p2 = p1;

%% Get the resulting positions
s_crit = determine_s_crit_1layer(kc_init1, w_init, z_crit, L_init, 0, 0.5, B, Binv);

positions = cell(size(lengths));
for i = 1:length(positions)
    positions{i} = position_kcp_2layer(B, Binv, kc_init1, kc_init2, 0, L_init, lengths(i), s_crit, p1, p2);
    
end

%% Plot the positions (3D)
f3d = figure('units','normalized','outerposition',[0 0 1 1]);
for i = length(positions):-1:1
    pi = positions{i};
    plot3(pi(3,:),pi(1,:),pi(2,:),'Displayname',sprintf("L = %d mm", lengths(i))); hold on;

end


% find the four corners of the tissue patch
[z_crit_found, idx] = min(abs(pi(3,:)-z_crit)); % find the z_crit index
S1 = [z_crit,        z_crit;
      pi(1,idx) - 5, pi(1,idx) + 5;
      pi(1,idx) - 5 , pi(2,idx) + 1];
S2 = S1;
S2(2,:) = S1(2,end:-1:1);
S = [S1(:,1) S2(:,1) S1(:,2) S2(:,2)];

patch(S(1,:), S(2,:), S(3,:), 'r','DisplayName', 'Tissue Boundary');

hold off;
title("Position Plots: 3D")
if length(lengths) < 7
    title("Position Plots: 2D")
    legend()
else
    title(sprintf("Position Plots (%d - %d mm | %d mm incs.): 3D", min(lengths), max(lengths), dL));
end
    
xlabel('z'); ylabel('x'); zlabel('y');
grid on;
axis equal;
view([-13.9832, 41.0820])

%% Plot the positions (2D)
f2d = figure('units','normalized','outerposition',[0 0 1 1]);
for i = length(positions):-1:1
    pi = positions{i};
    plot(pi(3,:),pi(2,:),'Displayname',sprintf("L = %d mm", lengths(i))); hold on;   
    
end

% plot the tissue boundary
xline(z_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');

hold off;

if length(lengths) < 7
    title("Position Plots: 2D")
    legend()
else
    title(sprintf("Position Plots (%d - %d mm | %d mm incs.): 2D", min(lengths), max(lengths), dL));
end

% annotate the layers
if abs(kc_init1) > abs(kc_init2)
    annotation('textbox',[0.2, 0.025, .2, .2], 'String','Hard Layer','FitBoxToText','On')
    annotation('textbox',[0.5, 0.025, .2, .2], 'String','Soft Layer','FitBoxToText','On')
    
elseif abs(kc_init1) < abs(kc_init2)
    annotation('textbox',[0.2, 0.025, .2, .2], 'String','Soft Layer','FitBoxToText','On')
    annotation('textbox',[0.5, 0.025, .2, .2], 'String','Hard Layer','FitBoxToText','On')
    
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
xline(z_crit, 'r--', 'Tissue Boundary');

if length(lengths) < 7
    title(sprintf("Mean Error | Overall mean: %.4f",...
        mean(norm_mean_error)));
    
else
    title(sprintf("Mean Error (%d - %d mm | %d mm incs.) | Overall mean: %.4f",...
        min(lengths), max(lengths), dL, mean(norm_mean_error)));
end

% annotate the layers
if abs(kc_init1) > abs(kc_init2)
    annotation('textbox',[0.2, 0.025, .2, .2], 'String','Hard Layer','FitBoxToText','On')
    annotation('textbox',[0.5, 0.025, .2, .2], 'String','Soft Layer','FitBoxToText','On')
    
elseif abs(kc_init1) < abs(kc_init2)
    annotation('textbox',[0.2, 0.025, .2, .2], 'String','Soft Layer','FitBoxToText','On')
    annotation('textbox',[0.5, 0.025, .2, .2], 'String','Hard Layer','FitBoxToText','On')
    
end

xlabel('z'); ylabel('Error (mm)');
grid on;
% ylim([0, 1]);
ylim([0,max([1, norm_mean_error])]);
% axis equal

%% save the plots
if ~save
    return
end

filename_pval = sprintf("kc-analysis_2layers_%dlengths_kc1_%.6f_kc2_%.6f_p1_%.6f_p2_%.6f",...
    length(lengths),kc_init1, kc_init2,  p1, p2);

% save the figure plots
savefig(f2d, directory + filename_pval + "_2d.fig")
disp("Saved figure: "+ directory + filename_pval + "_2d.fig");

savefig(f3d, directory + filename_pval + "_3d.fig");
disp("Saved figure: "+ directory + filename_pval + "_3d.fig");

savefig(ferr, directory + filename_pval + "_err.fig");
disp("Saved figure: "+ directory + filename_pval + "_err.fig");

% save the figures as PNG images.
saveas(f2d, directory + filename_pval + '_2d.png');
disp("Saved image: "+ directory + filename_pval + "_2d.png");

saveas(f3d, directory + filename_pval + '_3d.png');
disp("Saved image: "+ directory + filename_pval + "_3d.png");

saveas(ferr, directory + filename_pval + "_err.png");
disp("Saved image: "+ directory + filename_pval + "_err.png");
