%% ideal_singlebendinsertion_2020_0413.m
%
% used for the simulation of ideal single bend needle insertion.
%
% - written by Jin Seob Kim

clear all

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
ds = 0.5; % in mm
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
kc_s = 0.002;%0.0017; % for soft tissue (from IROS 2017 paper)
kc_h = 0.004; % for hard tissue (from IROS 2017 paper)

%% generate arclength array
ds = 0.5; % in mm
L1 = 60; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 90; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

%% optimal p value
%p = 0.55;
%p = 0.485;
p = 0.65; 

%% needle trajectories
kc1 = kc_s*(L/L1)^p;
w_init = [kc1;0;0];
[wv1,pmat1,Rmat1] = fn_intgEP_v1_1layer(w_init,kc1,0,0,ds,N1,B,Binv);
    
kc2 = kc_s*(L/L2)^p;
w_init = [kc2;0;0];
[wv2,pmat2,Rmat2] = fn_intgEP_v1_1layer(w_init,kc2,0,0,ds,N2,B,Binv);

kc3 = kc_s*(L/L3)^p;
w_init = [kc3;0;0];
[wv3,pmat3,Rmat3] = fn_intgEP_v1_1layer(w_init,kc3,0,0,ds,N3,B,Binv);

kc4 = kc_s*(L/L4)^p;
w_init = [kc4;0;0];
[wv4,pmat4,Rmat4] = fn_intgEP_v1_1layer(w_init,kc4,0,0,ds,N4,B,Binv);

kc5 = kc_s*(L/L5)^p;
w_init = [kc5;0;0];
[wv5,pmat5,Rmat5] = fn_intgEP_v1_1layer(w_init,kc5,0,0,ds,N5,B,Binv);

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
title(['p = ',num2str(p),', single-layer'])
legend('60mm','90mm','120mm','150mm','180mm')






    
