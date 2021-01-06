
%% main_doublebend_v2_tissue_90total.m
%
% old comment: use needleshape_fmincon_v2.m and data 09-02-2016
% double-bend, single homogeneous layer case
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

%% arclength (only consider full length insertion for now)
ds = 0.5; % in mm
L = 90;%80; % in mm
s = [0:ds:L];
N = length(s);

% location of FBG sensors (from the tip, in mm)
FBG3 = 11; % s_m3
FBG2 = 26; % s_m2
FBG1 = 70; % s_m1

sl1 = L - FBG1;
sl2 = L - FBG2;
sl3 = L - FBG3;

s_m1 = find(s == sl1);
s_m2 = find(s == sl2);
s_m3 = find(s == sl3);
s_m = [s_m1;s_m2;s_m3];

% arclength of turn (double bend)
s_turn = 40;
% s_shift = 5;
ix_s_turn = find(s == s_turn);% - s_shift);

s1 = s(1:ix_s_turn);
s2 = s(ix_s_turn:end);
N1 = length(s1);
N2 = length(s2);

%% load data
% New data on 04-25-2017 [1/m] (90 mm insertion, in tissue)
tot_num_exp = 10;

dbexpix = [2,4,6,8,10,12,14,16,18,20];

Data_raw_tot = zeros(6,20,tot_num_exp);
for i_exp = 1:tot_num_exp
    ii_exp = dbexpix(i_exp);
    sheet_num = ['T',num2str(ii_exp)];
    data_tmp = xlsread('changeframeResults_Insertion_0425_softphantom_steerat40mm_v1.xlsx',sheet_num,'C:V');
    Data_raw_tot(:,:,i_exp) = data_tmp(3:8,:);
end

% % data on 09-02-2016 [1/m] (90 mm insertion, in tissue)
% tot_num_exp = 10;
% 
% Data_raw_tot = zeros(6,20,tot_num_exp);
% for i_exp = 1:tot_num_exp
%     sheet_num = ['T',num2str(i_exp)];
%     data_tmp = xlsread('Results_insertion 0902.xlsx',sheet_num,'C:V');
%     Data_raw_tot(:,:,i_exp) = data_tmp(3:8,:);
% end

Data_raw_jig = squeeze(Data_raw_tot(:,end-1,:));
Data_raw_beam = squeeze(Data_raw_tot(:,end,:));

%% initialize
pmat_tot_jig = zeros(3,N,tot_num_exp);
pmat_tot_beam = zeros(3,N,tot_num_exp);

wmat_tot_jig = zeros(3,N,tot_num_exp);
wmat_tot_beam = zeros(3,N,tot_num_exp);

data_jig_tot = zeros(3,3,tot_num_exp);
data_beam_tot = zeros(3,3,tot_num_exp);

%% initial rotation angle (due to misalignment)
theta0 = 1.0;%0.7;
theta0 = theta0*pi/180;

%% total loop
kappa_c_mat = zeros(tot_num_exp,2);

for ii = 1:tot_num_exp

    disp(['%%%%%%%%%%%% number of experiment = ',num2str(ii),' %%%%%%%%%%%'])
    
    %% jig data
    disp(['%%%%% jig data %%%%%'])
    
    data_mat = Data_raw_jig(:,ii);
    w_m1 = -[data_mat(1);data_mat(2);0]; % equivalent to multiplying by Rot_z(-pi)
    w_m2 = -[data_mat(3);data_mat(4);0];
    w_m3 = -[data_mat(5);data_mat(6);0];
    data_jig = [w_m1,w_m2,w_m3]*1e-3; % in 1/mm

    data_jig_tot(:,:,ii) = data_jig;
    
    % parameters on intrinsic curvature (initial values)
    kc = data_jig(1,1)*(1+0.7);

    % initial configuration of a needle
    % initial guess of initial value eta
    eta = zeros(4,1);
    eta(1:3) = [kc;0;0];
    eta(4) = kc;

    % initial cost value
    scalef0 = 1;
    Cval = costfn_opt_v2_doublebend_v1(eta,data_jig,s_m,ix_s_turn,ds,N,B,Binv,scalef0);
    scalef = 1/Cval;

    % optimization
    tic; 
    maxiter = 1000;
    timelimit = 100;
    tolfun = 1e-6;%1e-4;

    oldopts = optimset('fmincon');
    %psoldopts = psoptimset;
    Tol = 1e-14*10^(ceil(log10(scalef)));

    options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
        'MaxFunEvals',10000,'Display','notify');
%     options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%         'MaxFunEvals',10000,'Display','iter','MaxIter',maxiter,'MaxFunEvals',100000);

    %psoptions = psoptimset(psoldopts,'Display','iter','MaxIter',maxiter);
    % gaset = saoptimset;
    % gaset = saoptimset(gaset,'Display','iter','TimeLimit',timelimit,'TolFun',tolfun);

    disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

    x0 = eta; % initial value
    LB = [-0.005*ones(3,1);0]; % lower bound
    UB = [0.005*ones(3,1);0.01]; % upper bound
    
    [x,fval,exitflag] = fmincon(@(x) costfn_opt_v2_doublebend_v1(x,data_jig,s_m,ix_s_turn,...
        ds,N,B,Binv,scalef),x0,[],[],[],[],LB,UB,[],options);

    t_f = toc;

%     x
%     fval
%     exitflag

    disp(['%>>>> optimization ends: time = ',num2str(t_f),' [sec]'])

    % final configuration
    kc_f = x(4);
    kappa_c_mat(ii,1) = kc_f;
    
    disp(['%>>>> kappa_c = ',num2str(kc_f)])

    disp(['%>> w_init = '])
    disp(num2str(x(1:3)))
    
    % intrinsic curvature kappa_0 (quadratic)
    kc_f1 = kc_f*(s1(end) - s1(1))/L;
    kc_f2 = kc_f*(s2(end) - s2(1))/L;
    
    k0_1_f = kc_f1*(s1(end)-s1).^2/L^2 - kc_f2*(1-s1(end)/L).*(1+s1(end)/L-2*s1/L);
    k0_2_f = -kc_f2*(1-s2/L).^2;
    k0_f = [k0_1_f(1:end-1),k0_2_f];

    k0prime1_f = -2*kc_f1/L*2*(s1(end) - s1) + kc_f2*(1-s1(end)/L).*(2/L);
    k0prime2_f = 2*kc_f2/L*(1 - s2/L);
    k0prime_f = [k0prime1_f(1:end-1),k0prime2_f];
    
    % intrinsic curvature \omega_0 and its derivative
    w0_f = [k0_f;zeros(size(s));zeros(size(s))];
    w0prime_f = [k0prime_f;zeros(size(s));zeros(size(s))];
    
%     w0_1_f = [k0_1_f;zeros(size(s1));zeros(size(s1))];
%     w0_2_f = [k0_2_f;zeros(size(s2));zeros(size(s2))];
% 
%     w0prime1_f = [k0prime1_f;zeros(size(s1));zeros(size(s1))];
%     w0prime2_f = [k0prime2_f;zeros(size(s2));zeros(size(s2))];
    
    % integration of the E-P equation
    wv_f = fn_intgEP_needle(x(1:3),w0_f,w0prime_f,ds,N,B,Binv);
    
%     wv_1_f = fn_intgEP_v1(x(1:3),w0_1_f,w0prime1_f,0,ds,N1,B,Binv);
%     wv_2_f = fn_intgEP_v1(-wv_1_f(:,end),w0_2_f,w0prime2_f,s2(1),ds,N2,B,Binv);
% 
%     wv_f = [wv_1_f(:,1:end-1),wv_2_f];
    
    % configuration
    Rmat = zeros(3,3,N);
    Rmat(:,:,1) = Rot_x(theta0);%eye(3);
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

    % save
    pmat_tot_jig(:,:,ii) = pmat;
    wmat_tot_jig(:,:,ii) = wv_f;
    
    %% beam data
    disp(['%%%%% beam data %%%%%'])

    data_mat = Data_raw_beam(:,ii);
    w_m1 = -[data_mat(1);data_mat(2);0]; % equivalent to multiplying by Rot_z(-pi)
    w_m2 = -[data_mat(3);data_mat(4);0];
    w_m3 = -[data_mat(5);data_mat(6);0];
    data_beam = [w_m1,w_m2,w_m3]*1e-3; % in 1/mm

    data_beam_tot(:,:,ii) = data_beam;

    % parameters on intrinsic curvature (initial values)
    kc = data_beam(1,1)*(1+0.7);

    % initial configuration of a needle
    % initial guess of initial value eta
    eta = zeros(4,1);
    eta(1:3) = [kc;0;0];
    eta(4) = kc;

    % initial cost value
    scalef0 = 1;
    Cval = costfn_opt_v2_doublebend_v1(eta,data_beam,s_m,ix_s_turn,ds,N,B,Binv,scalef0);
    scalef = 1/Cval;

    % optimization
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
    LB = [-0.005*ones(3,1);0]; % lower bound
    UB = [0.005*ones(3,1);0.01]; % upper bound
    
    [x,fval,exitflag] = fmincon(@(x) costfn_opt_v2_doublebend_v1(x,data_beam,s_m,ix_s_turn,...
        ds,N,B,Binv,scalef),x0,[],[],[],[],LB,UB,[],options);

    t_f = toc;
    
    % x
    % fval
    % exitflag

    disp(['%>>>> optimization ends: time = ',num2str(t_f),' [sec]'])

    % final configuration
    kc_f = x(4);
    kappa_c_mat(ii,2) = kc_f;
    
    disp(['%>>>> kappa_c = ',num2str(kc_f)])

    disp(['%>> w_init = '])
    disp(num2str(x(1:3)))
    
    % intrinsic curvature kappa_0 (quadratic)
    kc_f1 = kc_f*(s1(end) - s1(1))/L;
    kc_f2 = kc_f*(s2(end) - s2(1))/L;
    
    k0_1_f = kc_f1*(s1(end)-s1).^2/L^2 - kc_f2*(1-s1(end)/L).*(1+s1(end)/L-2*s1/L);
    k0_2_f = -kc_f2*(1-s2/L).^2;
    k0_f = [k0_1_f(1:end-1),k0_2_f];

    k0prime1_f = -2*kc_f1/L*2*(s1(end) - s1) + kc_f2*(1-s1(end)/L).*(2/L);
    k0prime2_f = 2*kc_f2/L*(1 - s2/L);
    k0prime_f = [k0prime1_f(1:end-1),k0prime2_f];
    
    % intrinsic curvature \omega_0 and its derivative
    w0_f = [k0_f;zeros(size(s));zeros(size(s))];
    w0prime_f = [k0prime_f;zeros(size(s));zeros(size(s))];
    
%     w0_1_f = [k0_1_f;zeros(size(s1));zeros(size(s1))];
%     w0_2_f = [k0_2_f;zeros(size(s2));zeros(size(s2))];
% 
%     w0prime1_f = [k0prime1_f;zeros(size(s1));zeros(size(s1))];
%     w0prime2_f = [k0prime2_f;zeros(size(s2));zeros(size(s2))];
    
    % integration of the E-P equation
    wv_f = fn_intgEP_needle(x(1:3),w0_f,w0prime_f,ds,N,B,Binv);
    
%     wv_1_f = fn_intgEP_v1(x(1:3),w0_1_f,w0prime1_f,0,ds,N1,B,Binv);
%     wv_2_f = fn_intgEP_v1(-wv_1_f(:,end),w0_2_f,w0prime2_f,s2(1),ds,N2,B,Binv);
% 
%     wv_f = [wv_1_f(:,1:end-1),wv_2_f];

    % configuration
    Rmat = zeros(3,3,N);
    Rmat(:,:,1) = Rot_x(theta0);%eye(3);
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

    % save
    pmat_tot_beam(:,:,ii) = pmat;
    wmat_tot_beam(:,:,ii) = wv_f;

    disp([' '])
end % end of total loop

%% extract needle trajectories from image analysis data for comparison
pos_img_tot = cell(tot_num_exp,1);
for jj = 1:tot_num_exp
    if jj <= 9
        sheet_num = ['0',num2str(jj)];
    elseif jj >= 10
        sheet_num = [num2str(jj)];
    end
    pos_img = xlsread('T0425_17.xlsx',sheet_num,'C:D');
    pos_img_tot{jj} = pos_img;
end
% last point is out-of-range

%% match image data with model results
pos_img_true_jig = cell(tot_num_exp,1);
pos_img_true_beam = cell(tot_num_exp,1);
for jj = 1:tot_num_exp
    if mod(jj,2) == 0
        pos_img_org = -pos_img_tot{jj};
    else
        pos_img_org = pos_img_tot{jj};
    end
    
    % select matched values
    ix_jig = find(pos_img_org(:,1) < pmat_tot_jig(3,end,jj));
    ix_beam = find(pos_img_org(:,1) < pmat_tot_beam(3,end,jj));
    
    pos_img_true_jig{jj} = pos_img_org(ix_jig,:);
    pos_img_true_beam{jj} = pos_img_org(ix_beam,:);
end

%% match image data with model results
pos_img_true_jig = cell(tot_num_exp,1);
pos_img_true_beam = cell(tot_num_exp,1);
for jj = 1:tot_num_exp
    pos_img_org = pos_img_tot{jj};
    
    % select matched values
    ix_jig = find(pos_img_org(:,1) < pmat_tot_jig(3,end,jj));
    ix_beam = find(pos_img_org(:,1) < pmat_tot_beam(3,end,jj));
    
    pos_img_true_jig{jj} = pos_img_org(ix_jig,:);
    pos_img_true_beam{jj} = pos_img_org(ix_beam,:);
end

%% tip deflection calculation
tip_defl_jig = zeros(tot_num_exp,1);
tip_defl_beam = zeros(tot_num_exp,1);

for jj = 1:tot_num_exp
    pos_img_t_jig = pos_img_true_jig{jj};
    pos_img_t_beam = pos_img_true_beam{jj};
    
    tip_defl1 = abs(pos_img_t_jig(end,2) - pmat_tot_jig(2,end,jj));
    tip_defl2 = abs(pos_img_t_beam(end,2) - pmat_tot_beam(2,end,jj));
    
    tip_defl_jig(jj) = tip_defl1;
    tip_defl_beam(jj) = tip_defl2;

    disp(['%%%%%%%%%%%% number of experiment = ',num2str(jj),' %%%%%%%%%%%'])
    disp(['difference of tip deflection = ',num2str(tip_defl1),' [mm] (jig)'])
    disp(['difference of tip deflection = ',num2str(tip_defl2),' [mm] (beam)'])
    disp([' '])
end

%% root mean square error of the trajectories
RMSE_jig = zeros(tot_num_exp,1);
RMSE_beam = zeros(tot_num_exp,1);

for jj = 1:tot_num_exp
    pos_img_t_jig = pos_img_true_jig{jj};
    pos_img_t_beam = pos_img_true_beam{jj};
    
    y_data_jig = pos_img_t_jig(:,2);
    z_data_jig = pos_img_t_jig(:,1);
    
    y_data_beam = pos_img_t_beam(:,2);
    z_data_beam = pos_img_t_beam(:,1);    
    
    pos_jig = pmat_tot_jig(:,:,jj);
    pos_beam = pmat_tot_beam(:,:,jj);
    
    RMSE_jig(jj) = extract_data_pos(y_data_jig,z_data_jig,pos_jig);
    RMSE_beam(jj) = extract_data_pos(y_data_beam,z_data_beam,pos_beam);

    disp(['%%%%%%%%%%%% number of experiment = ',num2str(jj),' %%%%%%%%%%%'])
    disp(['RMSE = ',num2str(RMSE_jig(jj)),' [mm] (jig)'])
    disp(['RMSE = ',num2str(RMSE_beam(jj)),' [mm] (beam)'])
    disp([' '])
end

%% plot curvature along x (tot_num_exp = 10 case)
figure;
for i_plot = 1:tot_num_exp
    subplot(5,2,i_plot)
    plot(s,wmat_tot_jig(1,:,i_plot),'k-','LineWidth',2)
    hold on
    plot(s,wmat_tot_beam(1,:,i_plot),'b-','LineWidth',2)
    plot(s(s_m),data_jig_tot(1,:,i_plot),'ko','MarkerSize',10,'LineWidth',2)
    plot(s(s_m),data_beam_tot(1,:,i_plot),'bs','MarkerSize',10,'LineWidth',2)
    hold off
    grid on
    if i_plot == 5 || i_plot == 6
        ylabel('{\omega}_x [1/mm]')
    end
    if i_plot ==9 || i_plot == 10
        xlabel('Arclength [mm]')
    end
    if i_plot == 10
        legend('model (jig)','model (beam)','measured (jig)','measured (beam)')
    end
end

%% plot curvature along y (tot_num_exp = 10 case)
figure;
for i_plot = 1:tot_num_exp
    subplot(5,2,i_plot)
    plot(s,wmat_tot_jig(2,:,i_plot),'k-','LineWidth',2)
    hold on
    plot(s,wmat_tot_beam(2,:,i_plot),'b-','LineWidth',2)
    plot(s(s_m),data_jig_tot(2,:,i_plot),'ko','MarkerSize',10,'LineWidth',2)
    plot(s(s_m),data_beam_tot(2,:,i_plot),'bs','MarkerSize',10,'LineWidth',2)
    hold off
    grid on
    if i_plot == 5 || i_plot == 6
        ylabel('{\omega}_y [1/mm]')
    end
    if i_plot ==9 || i_plot == 10
        xlabel('Arclength [mm]')
    end
    if i_plot == 10
        legend('model (jig)','model (beam)','measured (jig)','measured (beam)')
    end
end

%% plot torsion (tot_num_exp = 10 case)
figure;
for i_plot = 1:tot_num_exp
    subplot(5,2,i_plot)
    plot(s,wmat_tot_jig(3,:,i_plot),'k-','LineWidth',2)
    hold on
    plot(s,wmat_tot_beam(3,:,i_plot),'b-','LineWidth',2)
    hold off
    grid on
    if i_plot == 5 || i_plot == 6
        ylabel('{\omega}_z [1/mm]')
    end
    if i_plot ==9 || i_plot == 10
        xlabel('Arclength [mm]')
    end
    if i_plot == 10
        legend('model (jig)','model (beam)')
    end
end


%% plot needle trajectories (tot_num_exp = 10 case)
% y-z plane
figure;
for i_plot = 1:tot_num_exp    
    subplot(5,2,i_plot)
    plot(squeeze(pmat_tot_jig(3,:,i_plot)),squeeze(pmat_tot_jig(2,:,i_plot)),'k-','LineWidth',2)
    hold on
    plot(squeeze(pmat_tot_beam(3,:,i_plot)),squeeze(pmat_tot_beam(2,:,i_plot)),'b-','LineWidth',2)
    plot(pos_img_true_jig{i_plot}(1:2:end,1),pos_img_true_jig{i_plot}(1:2:end,2),'ko','MarkerSize',8)
    plot(pos_img_true_beam{i_plot}(1:2:end,1),pos_img_true_beam{i_plot}(1:2:end,2),'ko','MarkerSize',8)
    hold off
    grid on
    axis equal
    axis([0 100 -10 10])
    if i_plot == 5 || i_plot == 6
        ylabel('Y [mm]')
    end
    if i_plot ==9 || i_plot == 10
        xlabel('Z [mm]')
    end
    if i_plot == 10
        legend('model (jig)','model (beam)','image')
    end
end

% x-z plane
figure;
for i_plot = 1:tot_num_exp    
    subplot(5,2,i_plot)
    plot(squeeze(pmat_tot_jig(3,:,i_plot)),squeeze(pmat_tot_jig(1,:,i_plot)),'k-','LineWidth',2)
    hold on
    plot(squeeze(pmat_tot_beam(3,:,i_plot)),squeeze(pmat_tot_beam(1,:,i_plot)),'b-','LineWidth',2)
    hold off
    grid on
    axis equal
    axis([0 90 -5 5])
    if i_plot == 5 || i_plot == 6
        ylabel('X [mm]')
    end
    if i_plot ==9 || i_plot == 10
        xlabel('Z [mm]')
    end
    if i_plot == 10
        legend('model (jig)','model (beam)')
    end
end




% %% 3D plot
% figure;
% for i_plot = 1:10
%     subplot(5,2,i_plot)
%     plot3(pmat_tot_jig(1,:,i_plot),pmat_tot_jig(2,:,i_plot),pmat_tot_jig(3,:,i_plot),'k-','LineWidth',2)
%     hold on
%     plot3(pmat_tot_beam(1,:,i_plot),pmat_tot_beam(2,:,i_plot),pmat_tot_beam(3,:,i_plot),'b-','LineWidth',2)
% %     xlabel('X')
% %     ylabel('Y')
% %     zlabel('Z')
%     grid on
%     axis equal
%     if i_plot == 10
%         legend('model (jig)','model (beam)')
%     end
% end
% 
% figure;
% plot3(pmat_tot_jig(1,:,i_plot),pmat_tot_jig(2,:,i_plot),pmat_tot_jig(3,:,i_plot),'k-','LineWidth',2)
% grid on
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% axis equal
% axis([-5 0 -10 0 0 100])


