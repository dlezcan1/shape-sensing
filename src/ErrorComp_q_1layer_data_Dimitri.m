%% ErrorComp_q_1layer_data_Dimitri.m
%
% from ideal_singlebendinsertion_2020_0413.m
%
% used for the simulation of data-based single bend needle insertion.
%
% - written by Jin Seob Kim
% - edited  by Dimitri Lezcano

clear all;

% saving set-up
save_bool = true;
directory = "Dimitri/Data/";
directory = directory + "Archive/1-Layer/TipAError_Cost_Data/";
directory = directory + "kc_p_w_init_q/";
directory = directory + "exp_08/";
file_save = directory + "kc_singlelayer_p-%.3f_q-optim_data";

file_save = file_save + "_jig";
disp(directory)

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

%% intrinsic curvature
% first reference value when L = 90 mm
kc_s = 0.0025166;

% data
w_init_ref = [ 0.0036446; 0.00062854; -0.0099304 ];
theta0 = 1.0*pi/180;

%% generate arclength array
ds = 0.5; % in mm
L1 = 90; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 105; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

lengths = [L1 L2 L3 L4 L5];

%% q values as an array
p = 0.5860;
q_start = 0.0;
q_end = 1.0;
Nq = 101;
dp = (q_end-q_start)/(Nq-1);
qv = [q_start:dp:q_end];

%=========================================================================
%% main loop
total_error = zeros(1,Nq);
for i = 1:Nq
    q = qv(i);
    
    % needle trajectories
    kc1 = kc_s*(L/L1)^p;
    % w_init = [kc1;0;0];
    w_init = w_init_ref * (L/L1)^q;
    [wv1,pmat1,Rmat1] = fn_intgEP_v1_1layer(w_init,kc1,theta0,0,ds,N1,B,Binv);

    kc2 = kc_s*(L/L2)^p;
    % w_init = [kc2;0;0];
    w_init = w_init_ref * (L/L2)^q;
    [wv2,pmat2,Rmat2] = fn_intgEP_v1_1layer(w_init,kc2,theta0,0,ds,N2,B,Binv);

    kc3 = kc_s*(L/L3)^p;
    % w_init = [kc3;0;0];
    w_init = w_init_ref * (L/L3)^q;
    [wv3,pmat3,Rmat3] = fn_intgEP_v1_1layer(w_init,kc3,theta0,0,ds,N3,B,Binv);

    kc4 = kc_s*(L/L4)^p;
    % w_init = [kc4;0;0];
    w_init = w_init_ref * (L/L4)^q;
    [wv4,pmat4,Rmat4] = fn_intgEP_v1_1layer(w_init,kc4,theta0,0,ds,N4,B,Binv);

    kc5 = kc_s*(L/L5)^p;
    % w_init = [kc5;0;0];
    w_init = w_init_ref * (L/L5)^q;
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
q_best = qv(min_idx);
fprintf('Best q-value: %.5f\n', q_best)
disp(" ")

%% Generate q_best insertions
% needle trajectories
kc1 = kc_s*(L/L1)^p;
% w_init = [kc1;0;0];
w_init = w_init_ref * (L/L1)^q_best;
[wv1,pmat1,Rmat1] = fn_intgEP_v1_1layer(w_init,kc1,theta0,0,ds,N1,B,Binv);

kc2 = kc_s*(L/L2)^p;
% w_init = [kc2;0;0];
w_init = w_init_ref * (L/L2)^q_best;
[wv2,pmat2,Rmat2] = fn_intgEP_v1_1layer(w_init,kc2,theta0,0,ds,N2,B,Binv);

kc3 = kc_s*(L/L3)^p;
% w_init = [kc3;0;0];
w_init = w_init_ref * (L/L3)^q_best;
[wv3,pmat3,Rmat3] = fn_intgEP_v1_1layer(w_init,kc3,theta0,0,ds,N3,B,Binv);

kc4 = kc_s*(L/L4)^p;
% w_init = [kc4;0;0];
w_init = w_init_ref * (L/L4)^q_best;
[wv4,pmat4,Rmat4] = fn_intgEP_v1_1layer(w_init,kc4,theta0,0,ds,N4,B,Binv);

kc5 = kc_s*(L/L5)^p;
% w_init = [kc5;0;0];
w_init = w_init_ref * (L/L5)^q_best;
[wv5,pmat5,Rmat5] = fn_intgEP_v1_1layer(w_init,kc5,theta0,0,ds,N5,B,Binv);

% error
pmat_total = cell(1,5);
pmat_total{1} = pmat1;
pmat_total{2} = pmat2;
pmat_total{3} = pmat3;
pmat_total{4} = pmat4;
pmat_total{5} = pmat5;
    
%% plots
f1 = figure(1);
plot(qv,total_error,'ko-','linewidth',2)
xlabel('q')
ylabel('error [mm^2]')
title(sprintf("single-layer insertion: best p' = %.5f",q_best))
grid on

f2 = figure(2);
set(gcf, 'units','normalized','position',[0 0 1 1]);
f3 = figure(3);
set(gcf, 'units','normalized','position',[0 0 1 1]);
for i = length(pmat_total):-1:1
    figure(2);
    pmati = pmat_total{i};
    plot3(pmati(3,:), pmati(1,:), pmati(2,:), 'linewidth', 2, 'displayname', sprintf('L = %.1f', lengths(i)));
    hold on;
    
    figure(3);
    subplot(2,1,1);
    plot(pmati(3,:), pmati(1,:), 'linewidth', 2, 'displayname', sprintf('L = %.1f', lengths(i))); 
    hold on;
    
    subplot(2,1,2);
    plot(pmati(3,:), pmati(2,:), 'linewidth', 2, 'displayname', sprintf('L = %.1f', lengths(i))); 
    hold on;
    
end
figure(2);
hold off;
xlabel('z [mm]'); ylabel('x [mm]'); zlabel('y [mm]')
title(sprintf("3-D single-layer insertion: best p' = %.5f", q_best))
legend();
grid on; axis equal;

figure(3);
sgtitle(sprintf("2-D single-layer insertion: best p' = %.5f", q_best))
subplot(2,1,1);
hold off;
xlabel('z [mm]'); ylabel('x [mm]');
title('x-z axis view');
legend('location','northeastoutside');
grid on; axis equal;

subplot(2,1,2);
hold off;
xlabel('z [mm]'); ylabel('y [mm]');
title('y-z axis view');
grid on; axis equal;

%% saving
file_save = sprintf(file_save, p);

msg = "";
msg = msg + "fixed p-value and optimizing q-value\n";
msg = msg + sprintf('kc(90 mm), %.6f\n', kc_s);
msg = msg + sprintf('w_init(90 mm), [ %.5f; %.5f; %.5f ]\n', w_init_ref);
msg = msg + sprintf('p, %.5f\n', p);
msg = msg + sprintf('best q, %.5f\n', q_best);
msg = msg + sprintf('TipAError (mm^2), %.3f\n', min_err);

disp(msg);

if save_bool
    % figure
    savefig(f1, file_save + '_q-err.fig');
    fprintf('Saved figure: %s\n', file_save + '_q-err.fig');
    saveas(f1, file_save + '_q-err.png');
    fprintf('Saved image: %s\n', file_save + '_q-err.png');
    disp(" ");
    
    savefig(f2, file_save + '_3D.fig');
    fprintf('Saved figure: %s\n', file_save + '_3D.fig');
    saveas(f2, file_save + '_3D.png');
    fprintf('Saved image: %s\n', file_save + '_3D.png');
    disp(" ");
    
    savefig(f3, file_save + '_2D.fig');
    fprintf('Saved figure: %s\n', file_save + '_2D.fig');
    saveas(f3, file_save + '_2D.png');
    fprintf('Saved image: %s\n', file_save + '_2D.png');
    disp(" ");
    
    % csv file
    fid = fopen(file_save + '_err.csv', 'w');
    
    fprintf(fid, msg);
    
    fclose(fid);
    fprintf('Wrote file: %s\n', file_save + '_err.csv');
    disp(" ");
    
end   
