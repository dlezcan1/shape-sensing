%% inext_rod_singlebend_v3.m 
%
% solve inextensible rod equation for different insertion lengths
% use optimization to solve rod equation
% needs: fn_inextrod_w0.m, costfn_inextrod_w0.mm, cal_Rtraj.m
%
% inherited from inext_rod_sionglebend_v2.m --> conclusion:
% maintaining the same needle shape means that the forces
% applied to the needle should be decreasing.
% so I want to test reducing force values along the rod to have kappa_c
% relations.
%
% Conclusion: when applying p_p (order for load relation) about 2.5, then
% the resulting rod shapes look very similar. In that case kc ~ L^p where p
% is about -0.5. This is very close to what Dimitri obtained from his
% optmization work.
%
% - written by Jin Seob Kim

clear all

set(0,'DefaultAxesFontSize',28);

tic;

% e3 = [0;0;1];

%% preambles
% material properties
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

% arclength
ds = 0.5; % in mm
L = 90;%80; % in mm, reference insertion length
s = [0:ds:L];
N = length(s);

L1 = 60; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 90; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

%s_crit = 40;
%ix_s_crit = find(s == s_crit);

p_p = 2.5; % 2.5 is very close.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Case 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L1 = 60 mm
L1 = 60; N1 = L1/ds+1; s1 = linspace(0,L1,N1);

% force density values
press1 = 3e-2; %5e-3; %1e-2;%5e-3;
fric1 = 5e-3; %1e-3; %1e-3;%1e-3;%1e-3;%5e-3;%1e-3;

% external force
f1_ext = zeros(3,N1);
m1_ext = zeros(3,N1);

% single bend
f1_ext(2,:) = -press1;
f1_ext(3,:) = -fric1;

% rod equation solving
% initial values
% should satisfy s.t.
% \omega(0) = B^{-1} \int_0^L r x  R(s) f_ext ds
% f(0) = \int_0^L R(s) f_ext ds
% w(0) = [kc;0;0] and f(0) = [0;pr;fr]
% start with the results from the beam theory
kc = 0.05*press1*L1^2/2/BendStiff;
pr = 0.1*press1*L1;
fr = 0.1*fric1*L1;

x = [kc;pr;fr];

% initial cost value
scalef0 = 1;
Cval = costfn_inextrod_w0(x,m1_ext,f1_ext,B,Binv,ds,N1,scalef0);
scalef = 1/Cval;

% optimization
tic; 
maxiter = 500;
timelimit = 100;
tolfun = 1e-14;%1e-6;%1e-4;

oldopts = optimset('fmincon');
%psoldopts = psoptimset;
Tol = 1e-14*10^(ceil(log10(scalef)));

% options=optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',10000,'Display','notify');
options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',2000,'Display','iter','MaxIter',maxiter);

disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

x0 = x; % initial value

LB = [0;-50*pr;-50*fr]; % lower bound
UB = [50*kc;50*pr;50*fr]; % upper bound

[x,fval,exitflag] = fmincon(@(x) costfn_inextrod_w0(x,m1_ext,f1_ext,B,Binv,...
    ds,N1,scalef),x0,[],[],[],[],LB,UB,[],options);

t_f = toc;

% with friction
kc_f = x(1);
pr_f = x(2);
fr_f = x(3);
[wv1,fv1] = fn_inextrod_w0(kc_f,pr_f,fr_f,m1_ext,f1_ext,B,Binv,ds,N1);

% quadratic fit for single bend
y1 = wv1(1,1)/L1^2*(L1 - s1).^2; % --> very good even for large forces!

figure; 
plot(s1,wv1(1,:))
hold on
plot(s1,y1)
legend('\omega_1 (s)','\kappa_c (1 - s/L)^2')
grid on
xlabel('s [mm]')
ylabel('\omega_1')
title('L = 60 mm')
    
% configuration
[Rmat1,pmat1] = cal_Rtraj(wv1,ds,0);

% plot
% figure;
% plot(squeeze(pmat1(3,:)),squeeze(pmat1(2,:)),'k-','LineWidth',2)
% axis equal
% grid on
% xlabel('z [mm]')
% ylabel('y [mm]')

% figure;
% plot3(pmat(1,:),pmat(2,:),pmat(3,:),'k-','LineWidth',2)
% axis equal
% grid on
% xlabel('z [mm]')
% ylabel('y [mm]')
% zlabel('x [mm]')

% figure;
% plot(s1,wv1(1,:))
% %axis equal
% grid on
% xlabel('s [mm]')
% ylabel('\omega_x')
% 
% figure;
% plot(s1,fv(2,:))
% %axis equal
% grid on
% xlabel('s [mm]')
% ylabel('internal force')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L2 = 90 mm
L2 = 90; N2 = L2/ds+1; s2 = linspace(0,L2,N2);

% force density values
press2 = press1*(L1/L2)^p_p; %1e-2;%5e-3;
fric2 = fric1*(L1/L2)^p_p; % 2e-3; %1e-3;%1e-3;%1e-3;%5e-3;%1e-3;

% external force
f2_ext = zeros(3,N2);
m2_ext = zeros(3,N2);

% single bend
f2_ext(2,:) = -press2;
f2_ext(3,:) = -fric2;

% rod equation solving
% initial values
% should satisfy s.t.
% \omega(0) = B^{-1} \int_0^L r x  R(s) f_ext ds
% f(0) = \int_0^L R(s) f_ext ds
% w(0) = [kc;0;0] and f(0) = [0;pr;fr]
% start with the results from the beam theory
kc = 0.05*press2*L2^2/2/BendStiff;
pr = 0.1*press2*L2;
fr = 0.1*fric2*L2;

x = [kc;pr;fr];

% initial cost value
scalef0 = 1;
Cval = costfn_inextrod_w0(x,m2_ext,f2_ext,B,Binv,ds,N2,scalef0);
scalef = 1/Cval;

% optimization
tic; 
maxiter = 500;
timelimit = 100;
tolfun = 1e-14;%1e-6;%1e-4;

oldopts = optimset('fmincon');
%psoldopts = psoptimset;
Tol = 1e-14*10^(ceil(log10(scalef)));

% options=optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',10000,'Display','notify');
options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',2000,'Display','iter','MaxIter',maxiter);

disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

x0 = x; % initial value

LB = [0;-50*pr;-50*fr]; % lower bound
UB = [50*kc;50*pr;50*fr]; % upper bound

[x,fval,exitflag] = fmincon(@(x) costfn_inextrod_w0(x,m2_ext,f2_ext,B,Binv,...
    ds,N2,scalef),x0,[],[],[],[],LB,UB,[],options);

t_f = toc;

% with friction
kc_f = x(1);
pr_f = x(2);
fr_f = x(3);
[wv2,fv2] = fn_inextrod_w0(kc_f,pr_f,fr_f,m2_ext,f2_ext,B,Binv,ds,N2);

% quadratic fit for single bend
y2 = wv2(1,1)/L2^2*(L2 - s2).^2; % --> very good even for large forces!

figure; 
plot(s2,wv2(1,:))
hold on
plot(s2,y2)
legend('\omega_1 (s)','\kappa_c (1 - s/L)^2')
grid on
xlabel('s [mm]')
ylabel('\omega_1')
title('L = 90 mm')
    
% configuration
[Rmat2,pmat2] = cal_Rtraj(wv2,ds,0);

% % plot
% figure;
% plot(squeeze(pmat2(3,:)),squeeze(pmat2(2,:)),'k-','LineWidth',2)
% axis equal
% grid on
% xlabel('z [mm]')
% ylabel('y [mm]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Case 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L3 = 120 mm
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);

% force density values
press3 = press1*(L1/L3)^p_p; %1e-2;%5e-3;
fric3 = fric1*(L1/L3)^p_p; % 2e-3; %1e-3;%1e-3;%1e-3;%5e-3;%1e-3;

% external force
f3_ext = zeros(3,N3);
m3_ext = zeros(3,N3);

% single bend
f3_ext(2,:) = -press3;
f3_ext(3,:) = -fric3;

% rod equation solving
% initial values
% should satisfy s.t.
% \omega(0) = B^{-1} \int_0^L r x  R(s) f_ext ds
% f(0) = \int_0^L R(s) f_ext ds
% w(0) = [kc;0;0] and f(0) = [0;pr;fr]
% start with the results from the beam theory
kc = 0.05*press3*L3^2/2/BendStiff;
pr = 0.1*press3*L3;
fr = 0.1*fric3*L3;

x = [kc;pr;fr];

% initial cost value
scalef0 = 1;
Cval = costfn_inextrod_w0(x,m3_ext,f3_ext,B,Binv,ds,N3,scalef0);
scalef = 1/Cval;

% optimization
tic; 
maxiter = 500;
timelimit = 100;
tolfun = 1e-14;%1e-6;%1e-4;

oldopts = optimset('fmincon');
%psoldopts = psoptimset;
Tol = 1e-14*10^(ceil(log10(scalef)));

% options=optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',10000,'Display','notify');
options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',2000,'Display','iter','MaxIter',maxiter);

disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

x0 = x; % initial value

LB = [0;-50*pr;-50*fr]; % lower bound
UB = [50*kc;50*pr;50*fr]; % upper bound

[x,fval,exitflag] = fmincon(@(x) costfn_inextrod_w0(x,m3_ext,f3_ext,B,Binv,...
    ds,N3,scalef),x0,[],[],[],[],LB,UB,[],options);

t_f = toc;

% with friction
kc_f = x(1);
pr_f = x(2);
fr_f = x(3);
[wv3,fv3] = fn_inextrod_w0(kc_f,pr_f,fr_f,m3_ext,f3_ext,B,Binv,ds,N3);

% quadratic fit for single bend
y3 = wv3(1,1)/L3^2*(L3 - s3).^2; % --> very good even for large forces!

figure; 
plot(s3,wv3(1,:))
hold on
plot(s3,y3)
legend('\omega_1 (s)','\kappa_c (1 - s/L)^2')
grid on
xlabel('s [mm]')
ylabel('\omega_1')
title('L = 120 mm')
    
% configuration
[Rmat3,pmat3] = cal_Rtraj(wv3,ds,0);

% % plot
% figure;
% plot(squeeze(pmat3(3,:)),squeeze(pmat3(2,:)),'k-','LineWidth',2)
% axis equal
% grid on
% xlabel('z [mm]')
% ylabel('y [mm]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Case 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L4 = 150 mm
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);

% force density values
press4 = press1*(L1/L4)^p_p; %1e-2;%5e-3;
fric4 = fric1*(L1/L4)^p_p; % 2e-3; %1e-3;%1e-3;%1e-3;%5e-3;%1e-3;

% external force
f4_ext = zeros(3,N4);
m4_ext = zeros(3,N4);

% single bend
f4_ext(2,:) = -press4;
f4_ext(3,:) = -fric4;

% rod equation solving
% initial values
% should satisfy s.t.
% \omega(0) = B^{-1} \int_0^L r x  R(s) f_ext ds
% f(0) = \int_0^L R(s) f_ext ds
% w(0) = [kc;0;0] and f(0) = [0;pr;fr]
% start with the results from the beam theory
kc = 0.05*press4*L3^2/2/BendStiff;
pr = 0.1*press4*L3;
fr = 0.1*fric4*L3;

x = [kc;pr;fr];

% initial cost value
scalef0 = 1;
Cval = costfn_inextrod_w0(x,m4_ext,f4_ext,B,Binv,ds,N4,scalef0);
scalef = 1/Cval;

% optimization
tic; 
maxiter = 500;
timelimit = 100;
tolfun = 1e-14;%1e-6;%1e-4;

oldopts = optimset('fmincon');
%psoldopts = psoptimset;
Tol = 1e-14*10^(ceil(log10(scalef)));

% options=optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',10000,'Display','notify');
options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',2000,'Display','iter','MaxIter',maxiter);

disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

x0 = x; % initial value

LB = [0;-50*pr;-50*fr]; % lower bound
UB = [50*kc;50*pr;50*fr]; % upper bound

[x,fval,exitflag] = fmincon(@(x) costfn_inextrod_w0(x,m4_ext,f4_ext,B,Binv,...
    ds,N4,scalef),x0,[],[],[],[],LB,UB,[],options);

t_f = toc;

% with friction
kc_f = x(1);
pr_f = x(2);
fr_f = x(3);
[wv4,fv4] = fn_inextrod_w0(kc_f,pr_f,fr_f,m4_ext,f4_ext,B,Binv,ds,N4);

% quadratic fit for single bend
y4 = wv4(1,1)/L4^2*(L4 - s4).^2; % --> very good even for large forces!

figure; 
plot(s4,wv4(1,:))
hold on
plot(s4,y4)
legend('\omega_1 (s)','\kappa_c (1 - s/L)^2')
grid on
xlabel('s [mm]')
ylabel('\omega_1')
title('L = 150 mm')
    
% configuration
[Rmat4,pmat4] = cal_Rtraj(wv4,ds,0);

% % plot
% figure;
% plot(squeeze(pmat4(3,:)),squeeze(pmat4(2,:)),'k-','LineWidth',2)
% axis equal
% grid on
% xlabel('z [mm]')
% ylabel('y [mm]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   now compare rod shapes and kappa_c values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare
Lv = [L1, L2, L3, L4];
Kcv = [wv1(1,1), wv2(1,1), wv3(1,1), wv4(1,1)]; 

po = polyfit(log(Lv),log(Kcv),1) % determine p value in kc

figure;
plot(Lv, Kcv, 'ko-')
xlabel('insertion length [mm]')
ylabel('\kappa_c [1/mm]')

figure;
plot(squeeze(pmat1(3,:)),squeeze(pmat1(2,:)),'LineWidth',2)
hold on
plot(squeeze(pmat2(3,:)),squeeze(pmat2(2,:)),'LineWidth',2)
plot(squeeze(pmat3(3,:)),squeeze(pmat3(2,:)),'LineWidth',2)
plot(squeeze(pmat4(3,:)),squeeze(pmat4(2,:)),'LineWidth',2)
hold off
axis equal
grid on
xlabel('z [mm]')
ylabel('y [mm]')
legend('L = 60','L = 90','L = 120','L = 150')

