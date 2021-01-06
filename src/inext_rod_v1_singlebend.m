% 
% solve inextensible rod equation 
% use optimization
% needs: fn_inextrod_w0.m and costfn_inextrod_w0.m
%
% The results show that the assumed form of the intrinsic curvature is in
% an excellent match with simulated rod deformation case.
%
% - written by Jin Seob (Jesse) Kim

clear all

set(0,'DefaultAxesFontSize',28);

tic;

% e3 = [0;0;1];

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
ds = 0.5; % in mm
L = 90;%80; % in mm
s = [0:ds:L];
N = length(s);

%s_crit = 40;
%ix_s_crit = find(s == s_crit);

%% force density values
press = 1e-2;%5e-3;
fric = 1e-3;%1e-3;%1e-3;%5e-3;%1e-3;

%% external force
f_ext = zeros(3,N);
m_ext = zeros(3,N);

% single bend
f_ext(2,:) = -press;
f_ext(3,:) = -fric;

%% rod equation solving
% initial values
% should satisfy s.t.
% \omega(0) = B^{-1} \int_0^L r x  R(s) f_ext ds
% f(0) = \int_0^L R(s) f_ext ds
% w(0) = [kc;0;0] and f(0) = [0;pr;fr]
% start with the results from the beam theory
kc = 0.05*press*L^2/2/BendStiff;
pr = 0.1*press*L;
fr = 0.1*fric*L;

x = [kc;pr;fr];

% initial cost value
scalef0 = 1;
Cval = costfn_inextrod_w0(x,m_ext,f_ext,B,Binv,ds,N,scalef0);
scalef = 1/Cval;

% optimization
tic; 
maxiter = 500;
timelimit = 100;
tolfun = 1e-14;%1e-6;%1e-4;

oldopts = optimset('fmincon');
%psoldopts = psoptimset;
Tol = 1e-14*10^(ceil(log10(scalef)));

% options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',10000,'Display','notify');
options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',2000,'Display','iter','MaxIter',maxiter);

disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

x0 = x; % initial value

LB = [0;-50*pr;-50*fr]; % lower bound
UB = [50*kc;50*pr;50*fr]; % upper bound

[x,fval,exitflag] = fmincon(@(x) costfn_inextrod_w0(x,m_ext,f_ext,B,Binv,ds,N,scalef),...
    x0,[],[],[],[],LB,UB,[],options);

t_f = toc;

% with friction
kc_f = x(1);
pr_f = x(2);
fr_f = x(3);
[wv,fv] = fn_inextrod_w0(kc_f,pr_f,fr_f,m_ext,f_ext,B,Binv,ds,N);

% quadratic fit for single bend
y = wv(1,1)/L^2*(L - s).^2; % --> very good even for large forces!

figure; 
plot(s,wv(1,:))
hold on
plot(s,y)
legend('\omega_1 (s)','\kappa_c (1 - s/L)^2')
grid on
xlabel('s [mm]')
ylabel('\omega_1')
    
%% configuration
[Rmat,pmat] = cal_Rtraj(wv,ds,0);

%% plot
figure;
plot(squeeze(pmat(3,:)),squeeze(pmat(2,:)),'k-','LineWidth',2)
axis equal
grid on
xlabel('z [mm]')
ylabel('y [mm]')

figure;
plot3(pmat(1,:),pmat(2,:),pmat(3,:),'k-','LineWidth',2)
axis equal
grid on

figure;
plot(s,wv(1,:))
%axis equal
grid on

figure;
plot(s,fv(2,:))
%axis equal
grid on