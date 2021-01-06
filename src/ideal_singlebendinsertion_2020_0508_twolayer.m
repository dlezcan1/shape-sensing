%% ideal_singlebendinsertion_2020_0413_twolayer.m
%
% used for the simulation of ideal single bend needle insertion into two-layers.
%
% - written by Jin Seob Kim

clear all

global p % for simulation regarding p

set(0,'DefaultAxesFontSize',28);

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
ds = 0.5; % in mm, 1 mm seems fine as well
L = 90; % in mm
s = [0:ds:L];
N = length(s);

% % data points
% s_m1 = find(s == 20);
% s_m2 = find(s == 66);
% s_m3 = find(s == 79);
% s_m = [s_m1;s_m2;s_m3];

%% intrinsic curvature
% first reference value when L = 90 mm
kc_s = 0.002;%0.0017; % for soft tissue 
kc_h = 0.004; % for hard tissue

% % use this for hard-soft transition
% kc_s = 0.004;%0.0017; % for first tissue 
% kc_h = 0.002; % for second tissue

%% generate arclength array
L1 = 60; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 90; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

% z value for the boundary
z_crit = 45; % in mm

%% optimal p value
%p = 0.55;
%p = 0.49; % from Dimitri, optimization
p = 0.65; % from elastic rod equation solving (p = 0.65 is good)
%p = 0.795;

% p1 = p;%0.8;
% p2 = p;%0.658;

% a = -2/2;
% b = -1/2;

%% determine s_crit (from 90 mm homogeneous single bend insertion)
[w1,pos1] = fn_intgEP_v1_1layer([kc_s;0;0],kc_s,0,0,ds,N,B,Binv);
ix_crit_v = find(abs(pos1(3,:) - z_crit) <= ds/2); 
ix_crit = ix_crit_v(1);
s_crit = s(ix_crit);

%% needle trajectories
% kc_s_m = kc_s*(L/s_crit)^p*(s_crit/L1)^a;
% kc_h_m = kc_h*(L/(L-s_crit))^p*(1-(s_crit/L1)^2)^b;

kc_s_m = kc_s*(L/L1)^p;%(L/s_crit)^p1;%*(s_crit/L1)^a;
kc_h_m = kc_h*(L/L1)^p;%(L/(L-s_crit))^p2;%*(1-(s_crit/L1)^2)^b;
kc1 = kc_s_m*(s_crit/L1)^2 + kc_h_m*(1 - s_crit^2/L1^2);
w_init = [kc1;0;0];
[wv1,pmat1,Rmat1] = fn_intgEP_v2_2layers(w_init,kc_s,kc_h,z_crit,0,0,ds,N1,B,Binv);

kc_s_m = kc_s*(L/L2)^p;%(L/s_crit)^p1;%*(s_crit/L2)^a;
kc_h_m = kc_h*(L/L2)^p;%(L/(L-s_crit))^p2;%*(1-(s_crit/L2)^2)^b;
kc2 = kc_s_m*(s_crit/L2)^2 + kc_h_m*(1 - s_crit^2/L2^2);
w_init = [kc2;0;0];    
[wv2,pmat2,Rmat2] = fn_intgEP_v2_2layers(w_init,kc_s,kc_h,z_crit,0,0,ds,N2,B,Binv);

kc_s_m = kc_s*(L/L3)^p;%(L/s_crit)^p1;%*(s_crit/L3)^a;
kc_h_m = kc_h*(L/L3)^p;%(L/(L-s_crit))^p2;%*(1-(s_crit/L3)^2)^b;
kc3 = kc_s_m*(s_crit/L3)^2 + kc_h_m*(1 - s_crit^2/L3^2);
w_init = [kc3;0;0];
[wv3,pmat3,Rmat3] = fn_intgEP_v2_2layers(w_init,kc_s,kc_h,z_crit,0,0,ds,N3,B,Binv);

kc_s_m = kc_s*(L/L4)^p;%(L/s_crit)^p1;%*(s_crit/L4)^a;
kc_h_m = kc_h*(L/L4)^p;%(L/(L-s_crit))^p2;%*(1-(s_crit/L4)^2)^b;
kc4 = kc_s_m*(s_crit/L4)^2 + kc_h_m*(1 - s_crit^2/L4^2);
w_init = [kc4;0;0];
[wv4,pmat4,Rmat4] = fn_intgEP_v2_2layers(w_init,kc_s,kc_h,z_crit,0,0,ds,N4,B,Binv);

kc_s_m = kc_s*(L/L5)^p;%(L/s_crit)^p1;%*(s_crit/L5)^a;
kc_h_m = kc_h*(L/L5)^p;%(L/(L-s_crit))^p2;%*(1-(s_crit/L5)^2)^b;
kc5 = kc_s_m*(s_crit/L5)^2 + kc_h_m*(1 - s_crit^2/L5^2);
w_init = [kc5;0;0];
[wv5,pmat5,Rmat5] = fn_intgEP_v2_2layers(w_init,kc_s,kc_h,z_crit,0,0,ds,N5,B,Binv);

%% plots
figure;
plot(pmat1(3,:),pmat1(2,:))
hold on
plot(pmat2(3,:),pmat2(2,:))
plot(pmat3(3,:),pmat3(2,:))
plot(pmat4(3,:),pmat4(2,:))
plot(pmat5(3,:),pmat5(2,:))
axis equal
xlabel('z [mm]')
ylabel('y [mm]')
grid on
title(['p = ',num2str(p),', double-layer'])
legend('60mm','90mm','120mm','150mm','180mm')

%% error using area of a parameterized curve
x1 = pmat1(3,:); y1 = pmat1(2,:);
x2 = pmat2(3,:); y2 = pmat2(2,:);
x3 = pmat3(3,:); y3 = pmat3(2,:);
x4 = pmat4(3,:); y4 = pmat4(2,:);
x5 = pmat5(3,:); y5 = pmat5(2,:);

% indices for L's
ix60 = find(s2 == 60);
ix90 = find(s3 == 90);
ix120 = find(s4 == 120);
ix150 = find(s5 == 150);

% curve area error calculation
A_error = zeros(1,5);

% L = 60 case
A_error(1) = Acurve_param(-y2(1:ix60),x2(1:ix60)) - Acurve_param(-y1,x1) ...
    + Acurve_param(-y3(1:ix60),x3(1:ix60)) - Acurve_param(-y2(1:ix60),x2(1:ix60)) ...
    + Acurve_param(-y4(1:ix60),x4(1:ix60)) - Acurve_param(-y3(1:ix60),x3(1:ix60)) ...
    + Acurve_param(-y5(1:ix60),x5(1:ix60)) - Acurve_param(-y4(1:ix60),x4(1:ix60));
A_error(1) = A_error(1)/4;

% L = 90 case
A_error(2) = Acurve_param(-y3(1:ix90),x3(1:ix90)) - Acurve_param(-y2,x2) ...
    + Acurve_param(-y4(1:ix90),x4(1:ix90)) - Acurve_param(-y3(1:ix90),x3(1:ix90)) ...
    + Acurve_param(-y5(1:ix90),x5(1:ix90)) - Acurve_param(-y4(1:ix90),x4(1:ix90));
A_error(2) = A_error(2)/3;

% L = 120 case
A_error(3) = Acurve_param(-y4(1:ix120),x4(1:ix120)) - Acurve_param(-y3,x3) ...
    + Acurve_param(-y5(1:ix120),x5(1:ix120)) - Acurve_param(-y4(1:ix120),x4(1:ix120));
A_error(3) = A_error(2)/2;

% L = 150 case
A_error(4) = Acurve_param(-y5(1:ix150),x5(1:ix150)) - Acurve_param(-y4,x4);

% curve area error
A_er = abs(sum(A_error)/4);

% tip position error calculation
tip_error = zeros(1,5);

% L = 60 case
tip_error(1) = sqrt((x1(end) - x2(ix60))^2 + (y1(end) - y2(ix60))^2) ...
    + sqrt((x2(ix60) - x3(ix60))^2 + (y2(ix60) - y3(ix60))^2) ...
    + sqrt((x3(ix60) - x4(ix60))^2 + (y3(ix60) - y4(ix60))^2) ...
    + sqrt((x4(ix60) - x5(ix60))^2 + (y4(ix60) - y5(ix60))^2);
tip_error(1) = tip_error(1)/4*L1;

% L = 90 case
tip_error(2) = sqrt((x2(end) - x3(ix90))^2 + (y2(end) - y2(ix90))^2) ...
    + sqrt((x3(ix90) - x4(ix90))^2 + (y3(ix90) - y4(ix90))^2) ...
    + sqrt((x4(ix90) - x5(ix90))^2 + (y4(ix90) - y5(ix90))^2);
tip_error(2) = tip_error(2)/3*L2;

% L = 120 case
tip_error(3) = sqrt((x3(end) - x4(ix120))^2 + (y3(end) - y4(ix120))^2) ...
    + sqrt((x4(ix120) - x5(ix120))^2 + (y4(ix120) - y5(ix120))^2);
tip_error(3) = tip_error(3)/2*L3;

% L = 150 case
tip_error(4) = sqrt((x4(end) - x5(ix150))^2 + (y4(end) - y5(ix150))^2);
tip_error(4) = tip_error(4)*L4;

% total error
tip_er = sum(tip_error)/4;

disp(['>>>>> error = ',num2str(A_er + tip_er),', p = ',num2str(p)])

    
