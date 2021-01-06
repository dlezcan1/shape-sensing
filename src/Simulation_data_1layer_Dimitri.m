%% Simulation_data_1layer_Dimitri.m
%
% from ideal_singlebendinsertion_2020_0413.m
%
% used for the simulation of ideal single bend needle insertion.
%
% - written by Jin Seob Kim
% - edited by Dimitri Lezcano

clear all

%% material properties
% Stpainless Steel 304
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

%% arclength
% reference insertion (L = 90 mm)
ds = 0.5; % in mm
L = 90; % in mm

%% type of insertion
% data
w_init_ref = [ 0.0041315; -0.00073606; 0.0099075];
theta0 = 1.0*pi/180;

% ideal
ideal_insertion = false;

% p-value
p = 0.33;

%% intrinsic curvature
% first reference value when L = 90 mm
% kc_s = 0.002;%0.0017; % for soft tissue (from IROS 2017 paper)
%kc_h = 0.004; % for hard tissue (from IROS 2017 paper)
kc_s = 0.0030551;

%% File name for saving
directory = "Dimitri/Data";
directory = directory + "/Archive/1-Layer/TipAError_Cost_Data/";
file_name = directory + sprintf('kc_singlelayer_p_%.3f', p);

if ideal_insertion
    file_name = file_name + '_ideal';
    
else
    file_name = file_name + '_data';
    
end
    
%% generate arclength array
ds = 0.5; % in mm
L1 = 90; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 105; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

lengths = [L1, L2, L3, L4, L5];

%% Resulting Positions
if ideal_insertion
    w_init_ref = [kc_s; 0; 0];
    theta0 = 0;
end

% needle trajectories
kc1 = kc_s*(L/L1)^p;
% w_init = [kc1;0;0];
w_init = w_init_ref * (L/L1)^p;
[wv1,pmat1,Rmat1] = fn_intgEP_v1_1layer(w_init,kc1,theta0,0,ds,N1,B,Binv);

kc2 = kc_s*(L/L2)^p;
% w_init = [kc2;0;0];
w_init = w_init_ref * (L/L2)^p;
[wv2,pmat2,Rmat2] = fn_intgEP_v1_1layer(w_init,kc2,theta0,0,ds,N2,B,Binv);

kc3 = kc_s*(L/L3)^p;
% w_init = [kc3;0;0];
w_init = w_init_ref * (L/L3)^p;
[wv3,pmat3,Rmat3] = fn_intgEP_v1_1layer(w_init,kc3,theta0,0,ds,N3,B,Binv);

kc4 = kc_s*(L/L4)^p;
% w_init = [kc4;0;0];
w_init = w_init_ref * (L/L4)^p;
[wv4,pmat4,Rmat4] = fn_intgEP_v1_1layer(w_init,kc4,theta0,0,ds,N4,B,Binv);

kc5 = kc_s*(L/L5)^p;
% w_init = [kc5;0;0];
w_init = w_init_ref * (L/L5)^p;
[wv5,pmat5,Rmat5] = fn_intgEP_v1_1layer(w_init,kc5,theta0,0,ds,N5,B,Binv);

positions = {pmat1, pmat2, pmat3, pmat4, pmat5};

%% Plot of Position 3D
f3d = figure('units','normalized','outerposition',[0,0,1,1]);
for i = length(positions):-1:1
    pos = positions{i};
    
    plot3(pos(3,:), pos(1,:), pos(2,:), 'LineWidth',1.5, ...
        'DisplayName', sprintf("L = %d mm", lengths(i))); hold on;
    
end
title('3-D Plot of Trajectories');
xlabel('z [mm]'); ylabel('x [mm]'); zlabel('y [mm]');
lgd = legend();
lgd.FontSize = 10;
grid on;
axis equal;

if ~ideal_insertion
    savefig(f3d,file_name + "_3d.fig")
    disp("Saved figure: " + file_name + "_3d.fig");
    saveas(f3d,file_name + '_3d.png')
    disp("Saved figure: " + file_name + "_3d.png");
    disp(" ");
    
end

%% Plot of Positions 2D
f2d = figure('units','normalized','outerposition',[0,0,1,1]);
subplot(2,1,1);
for i = length(positions):-1:1
    pos = positions{i};
    
    plot(pos(3,:), pos(1,:), 'LineWidth',1.5, ...
        'DisplayName', sprintf("L = %d mm", lengths(i))); hold on;
    
end
title('z-x of Trajectories');
xlabel('z [mm]'); ylabel('x [mm]');
% legend()
grid on;
axis equal;
    
subplot(2,1,2);
for i = length(positions):-1:1
    pos = positions{i};
    
    plot(pos(3,:), pos(2,:), 'LineWidth',1.5, ...
        'DisplayName', sprintf("L = %d mm", lengths(i))); hold on;
    
end
title('z-y of Trajectories');
xlabel('z [mm]'); ylabel('y [mm]');
lgd = legend();
lgd.FontSize = 10;
grid on;
axis equal;

savefig(f2d,file_name + "_2d.fig")
disp("Saved figure: " + file_name + "_2d.fig");
saveas(f2d,file_name + '_2d.png')
disp("Saved figure: " + file_name + "_2d.png");
disp(" ");

%% Calculate Error

err = TipAerror(positions);
disp("Error Metric:");
disp(err);

%% Save Error
fid = fopen(file_name + '_err.csv', 'w');
fprintf(fid, 'p, %.4f\n', p);
fprintf(fid, 'w_init(90 mm), [ %.5f; %.5f; %.5f ]\n', w_init_ref);
fprintf(fid, 'TipAError (mm^2), %.5f\n', err);
fclose(fid);
fprintf('Wrote file: %s\n', file_name + '_err.csv');

