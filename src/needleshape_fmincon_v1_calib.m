%% needleshape_fmincon_v1_calib.m
%
% optimization with calibration data
%
% - written by Jin Seob (Jesse) Kim

clear all

set(0,'DefaultAxesFontSize',28);

tic;

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
s_m2 = find(s == 66);
s_m3 = find(s == 79);
s_m = [s_m1;s_m2;s_m3];

%% experimental data
data_mat = [...
    0	-0.053260596	-0.020481595	-0.009995091	-0.030276353	-0.01763456	-0.012049576;
0.5	0.510198657	0.208987204	0.116839581	-0.314382558	-0.204531504	-0.081774597;
1	1.108459172	0.771200149	0.513349288	-0.300154876	-0.08198115	-0.022067147;
1.4925	1.855114942	1.331274321	0.885596948	-0.272307571	-0.194855346	-0.077524051;
1.9802	2.268359175	1.723120163	1.085869891	-0.086945398	-0.117733612	-0.010100243;
2.4631	2.958731484	2.421058589	1.62479211	-0.112223155	-1.130563077	-0.427925648;
2.9412	3.111636385	2.79228902	1.594948492	-0.387014794	-0.301644554	-0.093470789];

data_mat = data_mat*1e-3;

data_cell = cell(size(data_mat,1),1);
for ii = 1:size(data_mat,1)
    w_m1 = [data_mat(ii,2);data_mat(ii,5);0];
    w_m2 = [data_mat(ii,3);data_mat(ii,6);0];
    w_m3 = [data_mat(ii,4);data_mat(ii,7);0];
    data_cell{ii} = [w_m1,w_m2,w_m3];
end

%% big loop to consider entire data
disp(['%%%%%%%%% start running the simulation %%%%%%%%%%'])
for ii = 1:size(data_cell,1)
    
%% parameters on intrinsic curvature (initial values)
kc = 0.0015;%abs(data_mat(ii,2)); 
    
%% initial configuration of a needle
% initial guess of initial value eta
eta = kc;
% eta = zeros(4,1);
% eta(1:3) = [kc;0;0];% + [0;0.001;0];
% eta(4) = kc;
% %eta(5) = alpha;

% angular deformation vector calculation
expno = ii; %input('number of experimental data = ');
data = data_cell{expno};

% single bend
k0 = kc*ones(size(s));%*(1 - alpha*s);
w0 = [k0;zeros(2,N)];
w0prime = zeros(3,N);

% k0prime = -kc*alpha*ones(1,N);
% w0prime = [k0prime;zeros(2,N)];

% w_init = eta(1:3);
% wv = fn_intgEP_v1(w_init,w0,w0prime,0,ds,N,B,Binv);
scalef0 = 1;
Cval = costfn_opt_v1_calib(eta,data,s_m,ds,N,B,Binv,scalef0);
scalef = 1/Cval;

%% optimization
maxiter = 1000;
timelimit = 100;
tolfun = 1e-4;

oldopts = optimset('fmincon');
%psoldopts = psoptimset;
Tol = 1e-12*10^(ceil(log10(scalef)));

options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
    'MaxFunEvals',100000);
%     'Display','iter','MaxIter',maxiter,'MaxFunEvals',100000);
%psoptions = psoptimset(psoldopts,'Display','iter','MaxIter',maxiter);
% gaset = saoptimset;
% gaset = saoptimset(gaset,'Display','iter','TimeLimit',timelimit,'TolFun',tolfun);

disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

x0 = eta;
LB = 0;
UB = 0.1;
% x0 = eta; % initial value
% LB = [-0.1*ones(3,1);0]; % lower bound
% UB = [0.1*ones(3,1);0.1]; % upper bound

[x,fval,exitflag] = fmincon(@(x) costfn_opt_v1_calib(x,data,s_m,ds,N,B,Binv,scalef),...
    x0,[],[],[],[],LB,UB,[],options);

% x
% fval
% exitflag

% disp(['%>>>> optimization ends: time = ',num2str(t_f),' [sec]'])

% %% final configuration
% kc_f = x(4);
% %alpha_f = x(5);
% 
% k0_f = kc_f*ones(size(s));%*(1 - alpha_f*s);
% w0_f = [k0_f;zeros(2,N)];
% w0prime_f = zeros(3,N);
% % k0prime_f = -kc_f*alpha_f*ones(1,N);
% % w0prime_f = [k0prime_f;zeros(2,N)];
% 
% wv_f = fn_intgEP_v1(x(1:3),w0_f,w0prime_f,0,ds,N,B,Binv);
% 
% % configuration
% Rmat = zeros(3,3,N);
% Rmat(:,:,1) = eye(3);
% pmat = zeros(3,N);
% for i = 2:N
%     % orientation
%     W = matr(1/2*(wv_f(:,i-1) + wv_f(:,i)));
%     Rmat(:,:,i) = Rmat(:,:,i-1)*expm(ds*W);
%     
%     % position
%     e3vec = squeeze(Rmat(:,3,1:i));
%     if i == 2
%         pmat(:,i) = pmat(:,i-1) + squeeze(Rmat(:,3,i))*ds;
%     else
%         pmat(:,i) = Simpson_vec_int(e3vec,ds);
%     end        
% end

%% compare constant curvature
curv_pred(ii,1) = data_mat(expno,1);
curv_cal(ii,1) = x;%x(4);
%error_curv = abs(curv_cal - curv_pred);

% disp(['%>>>> predefined intrinsic curvature = ',num2str(curv_pred)])
% disp(['%>>>> calculated intrinsic curvature = ',num2str(curv_cal)])
% disp(['%>>>> error = ',num2str(error_curv)])
% disp(['%%%%%%%%%%%%%%%%'])

%% end of big loop
end

t_f = toc;
disp(['%>>>> optimization ends: time = ',num2str(t_f),' [sec]'])
disp([' '])
disp(['predefined | calculated'])
disp(num2str([curv_pred,curv_cal]))

% %% plot the configuration
% figure;
% plot(pmat(3,:),pmat(2,:),'k-','LineWidth',2)
% xlabel('Z')
% ylabel('Y')
% grid on
% %axis([0 90 -5 5])
% % axis equal
% 
% figure; 
% plot3(pmat(1,:),pmat(2,:),pmat(3,:),'k-','LineWidth',2)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% grid on
% %axis([-5 5 -5 5 0 90])
% % axis equal
