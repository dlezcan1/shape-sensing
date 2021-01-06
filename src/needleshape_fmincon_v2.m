%% needleshape_fmincon_v2.m
%
% quadratic intrinsic curvature
%
% - written by Jin Seob (Jesse) Kim

clear all

set(0,'DefaultAxesFontSize',28);

tic;

%% select experiment type
Exp_No = 1; % 1: single bend, 2: double bend

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
L = 90; % in mm
s = [0:ds:L];
N = length(s);

% data points
s_m1 = find(s == 20);
s_m2 = find(s == 64); 
s_m3 = find(s == 79);
s_m = [s_m1;s_m2;s_m3];

%% experimental data
if Exp_No == 1
    % single bend
    data_mat = [...
    5.1224	0.001409425	0.000190242	-0.00012069	-6.56572E-05	-3.20545E-05	-3.28922E-05;
    5.1056	0.001399409	0.000214456	-8.58137E-05	6.71027E-05	7.09877E-05	1.86719E-05;
    5.276	0.001454885	0.000198112	-0.000114602	-0.000209859	2.802E-05	-1.11403E-05;
    5.3332	0.00149648	6.89087E-05	-0.000256192	4.06711E-05	8.22363E-05	-4.26215E-05;
    4.1023	0.001118877	9.74258E-05	-0.000281958	1.72744E-05	7.74086E-05	9.21884E-05;
    4.7046	0.001308601	8.55429E-05	-0.000379238	4.84661E-05	0.000196857	0.000386022;
    4.782	0.001319828	0.000166501	-0.000348302	-6.01577E-05	0.000185959	0.000305884];

    % data_mat = [...
    % 5.1224	0.001409425	0.000190242	0.00012069	-6.56572E-05	-3.20545E-05	-3.28922E-05;
    % 5.1056	0.001399409	0.000214456	8.58137E-05	6.71027E-05	7.09877E-05	1.86719E-05;
    % 5.276	0.001454885	0.000198112	0.000114602	-0.000209859	2.802E-05	-1.11403E-05;
    % 5.3332	0.00149648	6.89087E-05	0.000256192	4.06711E-05	8.22363E-05	-4.26215E-05;
    % 4.1023	0.001118877	9.74258E-05	0.000281958	1.72744E-05	7.74086E-05	9.21884E-05;
    % 4.7046	0.001308601	8.55429E-05	0.000379238	4.84661E-05	0.000196857	0.000386022;
    % 4.782	0.001319828	0.000166501	0.000348302	-6.01577E-05	0.000185959	0.000305884];

    data_cell = cell(size(data_mat,1),1);
    for ii = 1:size(data_mat,1)
        w_m1 = [data_mat(ii,2);data_mat(ii,5);0];
        w_m2 = [data_mat(ii,3);data_mat(ii,6);0];
        w_m3 = [data_mat(ii,4);data_mat(ii,7);0];
        data_cell{ii} = [w_m1,w_m2,w_m3];
    end

elseif Exp_No == 2 % #4 works well...
    % double bend
    data_mat = [...
    0.3927	40	2.91891E-07	0.000129137	-0.000128253	-0.000392237	6.26755E-05	6.76245E-05;
    0.15618	40	5.32303E-05	0.000251301	-7.25244E-05	-0.000193127	2.60953E-05	4.69712E-05;
    -0.46274	33	0.000248368	0.000209287	-9.2101E-05	-0.000335428	5.81734E-05	5.17259E-05;
    -1.4538	18	0.000578038	1.03718E-05	-0.000112264	-0.000477503	6.57657E-05	0.000182943;
    1.551	50	-0.000361853	0.000194346	-0.00011768	-0.000288384	3.53521E-05	4.7628E-05;
    1.2063	54	-0.000249729	0.000163054	-0.000178782	-6.69952E-05	0.000225356	0.000111511;
    -0.35476	35	0.000227138	0.000166777	-0.000216408	-0.000285765	5.23098E-06	-4.74507E-05;
    -1.2613	25	0.000500679	0.000184286	-0.000211056	-0.000530013	2.00051E-05	1.89279E-05;
    -0.78438	40	0.000359484	0.000173517	-0.000266081	-0.000369652	6.20405E-05	7.06025E-05;
    0.97686	49	-0.00013003	-0.000125402	-0.000407824	-9.57762E-05	0.000118697	5.4369E-05];    

    data_cell = cell(size(data_mat,1),1);
    for ii = 1:size(data_mat,1)
        w_m1 = [data_mat(ii,3);data_mat(ii,6);0];
        w_m2 = [data_mat(ii,4);data_mat(ii,7);0];
        w_m3 = [data_mat(ii,5);data_mat(ii,8);0];
        data_cell{ii} = [w_m1,w_m2,w_m3];
    end
end

%% parameters on intrinsic curvature (initial values)
kc = 0.001;%0.0015; 

%% initial configuration of a needle
% initial guess of initial value eta
eta = zeros(4,1);
eta(1:3) = [kc;0;0];% + [0;0.001;0];
eta(4) = kc;

% angular deformation vector calculation
disp(['%%%%%%%%%%%%%%%%%%%'])
expno = input('number of experimental data = ');
% expno = 1;
data = data_cell{expno};

if Exp_No == 1 % single bend
    % intrinsic curvature
    k0 = kc*(1 - s/L).^2;
    w0 = [k0;zeros(1,N);zeros(1,N)];

    k0prime = -2*kc/L*(1 - s/L);
    w0prime = [k0prime;zeros(1,N);zeros(1,N)];
    
    % initial cost value
    w_init = eta(1:3);
    wv = fn_intgEP_v1(w_init,w0,w0prime,0,ds,N,B,Binv);
    scalef0 = 1;
    Cval = costfn_opt_v2_singlebend(eta,data,s_m,ds,N,B,Binv,scalef0);
    scalef = 1/Cval;
elseif Exp_No == 2
    % double bend
    ix_s_turn = find(s == data_mat(expno,2));

    s1 = s(1:ix_s_turn-1);
    s2 = s(ix_s_turn:end);
    N1 = length(s1);
    N2 = length(s2);
    
    k0_1 = kc*(1-alpha*s1);
    k0_2 = -kc*(1-alpha*s2);
    
    w0_1 = [k0_1;zeros(size(s1));zeros(size(s1))];
    w0_2 = [k0_2;zeros(size(s2));zeros(size(s2))];
    
    k0prime1 = -kc*alpha*ones(size(s1));
    k0prime2 = kc*alpha*ones(size(s2));
    
    w0prime1 = [k0prime1;zeros(size(s1));zeros(size(s1))];
    w0prime2 = [k0prime2;zeros(size(s2));zeros(size(s2))];
    
    wv_1 = fn_intgEP_v1(eta(1:3),w0_1,w0prime1,0,ds,N1,B,Binv);
    wv_2 = fn_intgEP_v1(-wv_1(:,end),w0_2,w0prime2,s2(1),ds,N2,B,Binv);
    
    wv = [wv_1,wv_2];
%     k0 = zeros(1,N);
%     k0(1:ix_s_turn) = kc*(1 - alpha*s(1:ix_s_turn));
%     k0(ix_s_turn+1:end) = -kc*(1 - alpha*s(ix_s_turn+1:end));
%     w0 = [k0;zeros(2,N)];
% 
%     k0prime = zeros(1,N);
%     k0prime(1:ix_s_turn) = -kc*alpha;
%     k0prime(ix_s_turn+1:end) = kc*alpha;
%     w0prime = [k0prime;zeros(2,N)];
%     
%     w_init = eta(1:3);
%     wv = fn_intgEP_v1_doublebend(w_init,w0,w0prime,ix_s_turn,ds,N,B,Binv);
    scalef0 = 1;
    Cval = costfn_opt_v1_doublebend(eta,data,s_m,ix_s_turn,ds,N,B,Binv,scalef0);
    scalef = 1/Cval;
end

%% optimization
tic; 
maxiter = 1000;
timelimit = 100;
tolfun = 1e-6;%1e-4;

oldopts = optimset('fmincon');
%psoldopts = psoptimset;
Tol = 1e-14*10^(ceil(log10(scalef)));

options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',10000,'Display','notify');
%     'Display','iter','MaxIter',maxiter,'MaxFunEvals',100000);
    
%psoptions = psoptimset(psoldopts,'Display','iter','MaxIter',maxiter);
% gaset = saoptimset;
% gaset = saoptimset(gaset,'Display','iter','TimeLimit',timelimit,'TolFun',tolfun);

disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

x0 = eta; % initial value

limit_curv = 0.01;

LB = [-limit_curv*ones(3,1);0]; % lower bound
UB = [limit_curv*ones(3,1);limit_curv]; % upper bound

if Exp_No == 1
    [x,fval,exitflag] = fmincon(@(x) costfn_opt_v2_singlebend(x,data,s_m,ds,N,B,Binv,scalef),...
        x0,[],[],[],[],LB,UB,[],options);
elseif Exp_No == 2
    [x,fval,exitflag] = fmincon(@(x) costfn_opt_v2_doublebend(x,data,s_m,ix_s_turn,ds,N,B,Binv,scalef),...
        x0,[],[],[],[],LB,UB,[],options);
end

t_f = toc;

% x
% fval
% exitflag

disp(['%>>>> optimization ends: time = ',num2str(t_f),' [sec]'])

%% final configuration
kc_f = x(4);

disp(['%>>>> kappa_c = ',num2str(kc_f)])

if Exp_No == 1
    k0_f = kc_f*(1 - s/L).^2;
    w0_f = [k0_f;zeros(1,N);zeros(1,N)];

    k0prime_f = -2*kc_f/L*(1 - s/L);
    w0prime_f = [k0prime_f;zeros(1,N);zeros(1,N)];
    
    wv_f = fn_intgEP_v1(x(1:3),w0_f,w0prime_f,0,ds,N,B,Binv);
elseif Exp_No == 2
%     k0 = zeros(1,N);
%     k0(1:ix_s_turn) = kc_f*(1 - alpha_f*s(1:ix_s_turn));
%     k0(ix_s_turn+1:end) = -kc_f*(1 - alpha_f*s(ix_s_turn+1:end));
%     w0 = [k0;zeros(2,N)];
% 
%     k0prime = zeros(1,N);
%     k0prime(1:ix_s_turn) = -kc_f*alpha_f;
%     k0prime(ix_s_turn+1:end) = kc_f*alpha_f;
%     w0prime = [k0prime;zeros(2,N)];
    
    s1 = s(1:ix_s_turn-1);
    s2 = s(ix_s_turn:end);
    N1 = length(s1);
    N2 = length(s2);
    
    k0_1_f = kc_f*(1-alpha_f*s1);
    k0_2_f = -kc_f*(1-alpha_f*s2);
    
    w0_1_f = [k0_1_f;zeros(size(s1));zeros(size(s1))];
    w0_2_f = [k0_2_f;zeros(size(s2));zeros(size(s2))];
    
    k0prime1_f = -kc_f*alpha_f*ones(size(s1));
    k0prime2_f = kc_f*alpha_f*ones(size(s2));
    
    w0prime1_f = [k0prime1_f;zeros(size(s1));zeros(size(s1))];
    w0prime2_f = [k0prime2_f;zeros(size(s2));zeros(size(s2))];
    
    wv_f_1 = fn_intgEP_v1(x(1:3),w0_1_f,w0prime1_f,0,ds,N1,B,Binv);
    wv_f_2 = fn_intgEP_v1(-wv_f_1(:,end),w0_2_f,w0prime2_f,s2(1),ds,N2,B,Binv);
    
    wv_f = [wv_f_1,wv_f_2];
end

disp(['%>> w_init = '])
disp(num2str(x(1:3)))

% configuration
Rmat = zeros(3,3,N);
Rmat(:,:,1) = eye(3);
pmat = zeros(3,N);
for i = 2:N
    % orientation
    W = matr(1/2*(wv_f(:,i-1) + wv_f(:,i)));
    Rmat(:,:,i) = Rmat(:,:,i-1)*expm(ds*W);
    
    % position
    e3vec = squeeze(Rmat(:,3,1:i));
    if i == 2
        pmat(:,i) = pmat(:,i-1) + squeeze(Rmat(:,3,i))*ds;
    else
        pmat(:,i) = Simpson_vec_int(e3vec,ds);
    end        
end

%% compare tip deflection
y_tip_meas = abs(data_mat(expno,1));
y_tip_cal = abs(pmat(2,end));
error_tip = abs(y_tip_cal - y_tip_meas);
errorper = error_tip/y_tip_meas*100;

disp(['%>>>> measured tip deflection = ',num2str(y_tip_meas)])
disp(['%>>>> calculated tip deflection = ',num2str(y_tip_cal)])
disp(['%>>>> error = ',num2str(error_tip)])
%disp(['%>>>> error in % = ',num2str(errorper)])
disp(['%%%%%%%%%%%%%%%%'])

%% plot the configuration
figure; 
plot3(pmat(1,:),pmat(2,:),pmat(3,:),'k-','LineWidth',2)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
axis([-5 5 -5 5 0 90])
axis equal

figure;
subplot(3,1,1)
plot(pmat(3,:),pmat(2,:),'k-','LineWidth',2)
xlabel('Z [mm]')
ylabel('Y [mm]')
grid on
axis([0 90 -5 5])
axis equal

%% plot curvature
%figure;
subplot(3,1,2)
plot(s,wv_f(1,:),'k-','LineWidth',2)
hold on
plot(s(s_m),data(1,:),'ko','MarkerSize',10,'LineWidth',2)
hold off
grid on
xlabel('Arclength [mm]')
ylabel('{\omega}_x [1/mm]')
legend('model','measured')

%figure;
subplot(3,1,3)
plot(s,wv_f(2,:),'k-','LineWidth',2)
hold on
plot(s(s_m),data(2,:),'ko','MarkerSize',10,'LineWidth',2)
hold off
grid on
xlabel('Arclength [mm]')
ylabel('{\omega}_y [1/mm]')
legend('model','measured')

figure;
plot(s,wv_f(3,:),'k-','LineWidth',2)
grid on
xlabel('Arclength [mm]')
ylabel('{\omega}_z [1/mm]')