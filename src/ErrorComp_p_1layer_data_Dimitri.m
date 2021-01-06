%% ErrorComp_p_1layer_data_Dimitri.m
%
% from ideal_singlebendinsertion_2020_0413.m
%
% used for the simulation of data-based single bend needle insertion.
%
% - written by Jin Seob Kim
% - edited  by Dimitri Lezcano

clear all;

% saving set-up
save_bool = false;
directory = "Dimitri/Data/";
directory = directory + "Archive/1-Layer/TipAError_Cost_Data/";
directory = directory + "exp_01/";
file_save = directory + "kc_singlelayer_p-optim_data";

file_save = file_save + "_beam";

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

%% arclength
% reference insertion (L = 90 mm)
ds = 0.5; % in mm
L = 90; % in mm
% s = [0:ds:L];
% N = length(s);

% % data points
% s_m1 = find(s == 20);
% s_m2 = find(s == 66);
% s_m3 = find(s == 79);
% s_m = [s_m1;s_m2;s_m3];

%% intrinsic curvature
% first reference value when L = 90 mm
kc_s = 0.002;

% data
w_init_ref = [kc_s; 0; 0];
theta0 = 1.0*pi/180;

%% generate arclength array
ds = 0.5; % in mm
L1 = 90; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 105; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

%% p values as an array
p_start = 0.0;
p_end = 1.0;
Np = 101;
dp = (p_end-p_start)/(Np-1);
pv = [p_start:dp:p_end];

%=========================================================================
%% main loop
total_error = zeros(1,Np);
for i = 1:Np
    p = pv(i);
    
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
    % error
    pmat_total = cell(1,5);
    pmat_total{1} = pmat1;
    pmat_total{2} = pmat2;
    pmat_total{3} = pmat3;
    pmat_total{4} = pmat4;
    pmat_total{5} = pmat5;

    total_error(i) = TipAerror(pmat_total);
end

%% best p-value
[min_err, min_idx] = min(total_error);
p_best = pv(min_idx);
fprintf('Best p-value: %.5f\n', p_best)
disp(" ")

%% plots
f1 = figure();
plot(pv,total_error,'ko-','linewidth',2)
xlabel('p')
ylabel('error [mm^2]')
title(sprintf('single-layer insertion: best p = %.4f',p_best))
grid on

%% saving
if save_bool
    % figure
    savefig(f1, file_save + '_p-err.fig');
    fprintf('Saved figure: %s\n', file_save + '_p-err.fig');
    saveas(f1, file_save + '_p-err.png');
    fprintf('Saved image: %s\n', file_save + '_p-err.pmg');

    % csv file
    fid = fopen(file_save + '_err.csv', 'w');
    fprintf(fid, 'kc(90 mm), %.6f\n', kc_s);
    fprintf(fid, 'w_init(90 mm), [ %.5f; %.5f; %.5f ]\n', w_init_ref);
    fprintf(fid, 'best p, %.5f\n', p_best);
    fprintf(fid, 'TipAError (mm^2), %.3f\n', min_err);
    fclose(fid);
    fprintf('Wrote file: %s\n', file_save + '_err.csv');

else
    fprintf('kc(90 mm), %.6f\n', kc_s);
    fprintf('w_init(90 mm), [ %.5f; %.5f; %.5f ]\n', w_init_ref);
    fprintf('best p, %.5f\n', p_best);
    fprintf('TipAError (mm^2), %.3f\n', min_err);
    
end