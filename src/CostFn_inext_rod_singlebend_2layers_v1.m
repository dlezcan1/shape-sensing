function C = CostFn_inext_rod_singlebend_2layers_v1(qv,scalef_c)
%
% from inext_rod_singlebend_v4.m 
%
% solve inextensible rod equation for different insertion lengths
% use optimization to solve rod equation
% needs: fn_inextrod_w0.m, costfn_inextrod_w0.mm, cal_Rtraj.m
%
%
%
% - written by Jin Seob Kim

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

s_crit = 45;

% input
q1 = qv(1);
q2 = qv(2);

% array for kappa_c1 and kappa_c2
kc1vec = zeros(1,4);
kc2vec = zeros(1,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Case 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L1 = 60 mm
L1 = 60; N1 = L1/ds+1; s1 = linspace(0,L1,N1);

% force density values
press1_1 = 2e-2;%5e-3; %1e-2;%5e-3;
fric1_1 = 6e-3;%1e-3; %1e-3; %1e-3;%1e-3;%1e-3;%5e-3;%1e-3;
press1_2 = 1e-2;%1e-2;
fric1_2 = 3e-3;%2.5e-3;


% external force
f1_ext = zeros(3,N1);
m1_ext = zeros(3,N1);

ix_s1 = find(s1 == s_crit);

% single bend
f1_ext(2,1:ix_s1(1)) = -press1_1;
f1_ext(3,1:ix_s1(1)) = -fric1_1;
f1_ext(2,ix_s1(1)+1:end) = -press1_2;
f1_ext(3,ix_s1(1)+1:end) = -fric1_2;

% rod equation solving
% initial values
% should satisfy s.t.
% \omega(0) = B^{-1} \int_0^L r x  R(s) f_ext ds
% f(0) = \int_0^L R(s) f_ext ds
% w(0) = [kc;0;0] and f(0) = [0;pr;fr]
% start with the results from the beam theory
kc = 0.05*press1_1*L1^2/2/BendStiff;
pr = 0.1*press1_1*L1;
fr = 0.1*fric1_1*L1;

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

options=optimset(oldopts,'Algorithm','sqp','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',10000,'Display','off');
% options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',2000,'Display','iter','MaxIter',maxiter);

%disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

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

% extract kc1 and kc2 anbd compare profile
kc2vec(1) = wv1(1,ix_s1)/(1 - s_crit/L1)^2;
kc1vec(1) = (L1/s_crit)^2*(wv1(1,1) - kc2vec(1)*(1 - (s_crit/L1)^2));
y1 = zeros(1,N1);
for i = 1:ix_s1
    y1(i) = kc1vec(1)*((s_crit-s1(i))/L1)^2 + ...
        kc2vec(1)*(1-s_crit/L1)*(1+s_crit/L1-2*s1(i)/L1);
end
for i = ix_s1+1:N1
    y1(i) = kc2vec(1)*(1-s1(i)/L1)^2;
end
    
% configuration
[Rmat1,pmat1] = cal_Rtraj(wv1,ds,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L2 = 90 mm
L2 = 90; N2 = L2/ds+1; s2 = linspace(0,L2,N2);

% force density values
press2_1 = press1_1*(L1/L2)^q1; %1e-2;%5e-3;
fric2_1 = fric1_1*(L1/L2)^q1; % 2e-3; %1e-3;%1e-3;%1e-3;%5e-3;%1e-3;
press2_2 = press1_2*(L1/L2)^q2;
fric2_2 = fric1_2*(L1/L2)^q2;

% external force
f2_ext = zeros(3,N2);
m2_ext = zeros(3,N2);

ix_s2 = find(s2 == s_crit);

% single bend
f2_ext(2,1:ix_s2(1)) = -press2_1;
f2_ext(3,1:ix_s2(1)) = -fric2_1;
f2_ext(2,ix_s2(1)+1:end) = -press2_2;
f2_ext(3,ix_s2(1)+1:end) = -fric2_2;

% rod equation solving
% initial values
% should satisfy s.t.
% \omega(0) = B^{-1} \int_0^L r x  R(s) f_ext ds
% f(0) = \int_0^L R(s) f_ext ds
% w(0) = [kc;0;0] and f(0) = [0;pr;fr]
% start with the results from the beam theory
kc = 0.05*(press2_1 + press2_2)*L2^2/2/BendStiff;
pr = 0.1*(press2_1 + press2_2)*L2;
fr = 0.1*(fric2_1 + fric2_2)*L2;

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

options=optimset(oldopts,'Algorithm','sqp','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',10000,'Display','off');
% options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',2000,'Display','iter','MaxIter',maxiter);

%disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

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

% extract kc1 and kc2 anbd compare profile
kc2vec(2) = wv2(1,ix_s2)/(1 - s_crit/L2)^2;
kc1vec(2) = (L2/s_crit)^2*(wv2(1,1) - kc2vec(2)*(1 - (s_crit/L2)^2));
y2 = zeros(1,N2);
for i = 1:ix_s2
    y2(i) = kc1vec(2)*((s_crit-s2(i))/L2)^2 + ...
        kc2vec(2)*(1-s_crit/L2)*(1+s_crit/L2-2*s2(i)/L2);
end
for i = ix_s2+1:N2
    y2(i) = kc2vec(2)*(1-s2(i)/L2)^2;
end
    
% configuration
[Rmat2,pmat2] = cal_Rtraj(wv2,ds,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Case 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L3 = 120 mm
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);

% force density values
press3_1 = press1_1*(L1/L3)^q1; %1e-2;%5e-3;
fric3_1 = fric1_1*(L1/L3)^q1; % 2e-3; %1e-3;%1e-3;%1e-3;%5e-3;%1e-3;
press3_2 = press1_2*(L1/L3)^q2;
fric3_2 = fric1_2*(L1/L3)^q2;

% external force
f3_ext = zeros(3,N3);
m3_ext = zeros(3,N3);

ix_s3 = find(s3 == s_crit);

% single bend
f3_ext(2,1:ix_s3(1)) = -press3_1;
f3_ext(3,1:ix_s3(1)) = -fric3_1;
f3_ext(2,ix_s3(1)+1:end) = -press3_2;
f3_ext(3,ix_s3(1)+1:end) = -fric3_2;

% rod equation solving
% initial values
% should satisfy s.t.
% \omega(0) = B^{-1} \int_0^L r x  R(s) f_ext ds
% f(0) = \int_0^L R(s) f_ext ds
% w(0) = [kc;0;0] and f(0) = [0;pr;fr]
% start with the results from the beam theory
kc = 0.05*(press3_1 + press3_2)*L3^2/2/BendStiff;
pr = 0.1*(press3_1 + press3_2)*L3;
fr = 0.1*(fric3_1 + fric3_2)*L3;

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

options=optimset(oldopts,'Algorithm','sqp','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',10000,'Display','off');
% options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',2000,'Display','iter','MaxIter',maxiter);

%disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

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

% extract kc1 and kc2 anbd compare profile
kc2vec(3) = wv3(1,ix_s3)/(1 - s_crit/L3)^2;
kc1vec(3) = (L3/s_crit)^2*(wv3(1,1) - kc2vec(3)*(1 - (s_crit/L3)^2));
y3 = zeros(1,N3);
for i = 1:ix_s3
    y3(i) = kc1vec(3)*((s_crit-s3(i))/L3)^2 + ...
        kc2vec(3)*(1-s_crit/L3)*(1+s_crit/L3-2*s3(i)/L3);
end
for i = ix_s3+1:N3
    y3(i) = kc2vec(3)*(1-s3(i)/L3)^2;
end
    
% configuration
[Rmat3,pmat3] = cal_Rtraj(wv3,ds,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Case 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L4 = 150 mm
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);

% force density values
press4_1 = press1_1*(L1/L4)^q1; %1e-2;%5e-3;
fric4_1 = fric1_1*(L1/L4)^q1; % 2e-3; %1e-3;%1e-3;%1e-3;%5e-3;%1e-3;
press4_2 = press1_2*(L1/L4)^q2;
fric4_2 = fric1_2*(L1/L4)^q2;

% external force
f4_ext = zeros(3,N4);
m4_ext = zeros(3,N4);

ix_s4 = find(s4 == s_crit);

% single bend
f4_ext(2,1:ix_s4(1)) = -press4_1;
f4_ext(3,1:ix_s4(1)) = -fric4_1;
f4_ext(2,ix_s4(1)+1:end) = -press4_2;
f4_ext(3,ix_s4(1)+1:end) = -fric4_2;

% rod equation solving
% initial values
% should satisfy s.t.
% \omega(0) = B^{-1} \int_0^L r x  R(s) f_ext ds
% f(0) = \int_0^L R(s) f_ext ds
% w(0) = [kc;0;0] and f(0) = [0;pr;fr]
% start with the results from the beam theory
kc = 0.05*press4_1*L4^2/2/BendStiff;
pr = 0.1*press4_1*L4;
fr = 0.1*fric4_1*L4;

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

options=optimset(oldopts,'Algorithm','sqp','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',10000,'Display','off');
% options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',2000,'Display','iter','MaxIter',maxiter);

%disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

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

% extract kc1 and kc2 anbd compare profile
kc2vec(4) = wv4(1,ix_s4)/(1 - s_crit/L4)^2;
kc1vec(4) = (L4/s_crit)^2*(wv4(1,1) - kc2vec(4)*(1 - (s_crit/L4)^2));
y4 = zeros(1,N4);
for i = 1:ix_s4
    y4(i) = kc1vec(4)*((s_crit-s4(i))/L4)^2 + ...
        kc2vec(4)*(1-s_crit/L4)*(1+s_crit/L4-2*s4(i)/L4);
end
for i = ix_s4+1:N4
    y4(i) = kc2vec(4)*(1-s4(i)/L4)^2;
end

% configuration
[Rmat4,pmat4] = cal_Rtraj(wv4,ds,0);


%% =======================================================================
%                  Cost function value as the error
% ========================================================================
pmat_total = cell(1,4);
pmat_total{1} = pmat1;
pmat_total{2} = pmat2;
pmat_total{3} = pmat3;
pmat_total{4} = pmat4;

C = TipAerror(pmat_total)*scalef_c;

end
