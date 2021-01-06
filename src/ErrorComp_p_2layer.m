%% ErrorComp_p_2layer.m

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
kc_s = 0.002;%0.002;%0.0017; % for soft tissue 
kc_h = 0.004;%0.004; % for hard tissue

% % use this for hard-soft transition
kc_s = 0.004;%0.0017; % for first tissue 
kc_h = 0.002; % for second tissue

%% generate arclength array
L1 = 60; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 90; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

% z value for the boundary
z_crit = 35;%45; % in mm

%% determine s_crit (from 90 mm homogeneous single bend insertion)
[w1,pos1] = fn_intgEP_v1_1layer([kc_s;0;0],kc_s,0,0,ds,N,B,Binv);
ix_crit_v = find(abs(pos1(3,:) - z_crit) <= ds/2); 
ix_crit = ix_crit_v(1);
s_crit = s(ix_crit);

%% p values as an array
p_start = 0.4;
p_end = 0.7;
Np = 51;
dp = (p_end-p_start)/(Np-1);
pv = [p_start:dp:p_end];

%=========================================================================
%% main loop
total_error = zeros(1,Np);
for i = 1:Np
    p = pv(i);

    % needle trajectories
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
    
    % error
    pmat_total = cell(1,5);
    pmat_total{1} = pmat1;
    pmat_total{2} = pmat2;
    pmat_total{3} = pmat3;
    pmat_total{4} = pmat4;
    pmat_total{5} = pmat5;

    total_error(i) = TipAerror(pmat_total);
end

%% plots
figure;
plot(pv,total_error,'ko-','linewidth',2)
xlabel('p')
ylabel('error [mm^2]')
title('double-layer insertion')
grid on