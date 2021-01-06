%% main_singlebend_v2_2tissuelayers_90total_initangopt.m
%
% old comment: use needleshape_fmincon_v2.m and data 09-02-2016
% here kc are from individual homogeneous tissue insertion.
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

% z value for the boundary
z_crit = 45; % in mm

%% load data
% New data on 09-19-2016 two-layered phantom [1/m] (90 mm insertion, in tissue)
tot_num_exp = 10;

Data_raw_tot = zeros(6,20,tot_num_exp);
for i_exp = 1:tot_num_exp
    sheet_num = ['T',num2str(i_exp)];
%     data_tmp = xlsread('Results_Insertion 0920_TwoLayerPhantom_n180.xlsx',sheet_num,'C:V');    
    data_tmp = readmatrix('Results_Insertion 0920_TwoLayerPhantom_n180.xlsx','Sheet',sheet_num,'Range','C:V');
    Data_raw_tot(:,:,i_exp) = data_tmp(1:6,:);
end

Data_raw_jig = squeeze(Data_raw_tot(:,end-1,:));
Data_raw_beam = squeeze(Data_raw_tot(:,end,:));

%% initialize
wmat_tot_jig = zeros(3,N,tot_num_exp);
wmat_tot_beam = zeros(3,N,tot_num_exp);

data_jig_tot = zeros(3,3,tot_num_exp);
data_beam_tot = zeros(3,3,tot_num_exp);

%% kappa_c values from homogeneous tissue needle insertion (90 mm)
kc1_jig = 0.0017;
kc1_beam = 0.0020;

kc2_jig = 0.0027;
kc2_beam = 0.0032;

theta0 = 1.0*pi/180;
for ii = 1:tot_num_exp

    disp(['%%%%%%%%%%%% number of experiment = ',num2str(ii),' %%%%%%%%%%%'])
    
    %% jig data
    disp(['%%%%% jig data %%%%%'])
    
    data_mat = Data_raw_jig(:,ii);
    w_m1 = [data_mat(1);data_mat(2);0];
    w_m2 = [data_mat(3);data_mat(4);0];
    w_m3 = [data_mat(5);data_mat(6);0];
    data_jig = [w_m1,w_m2,w_m3]*1e-3; % in 1/mm

    data_jig_tot(:,:,ii) = data_jig;

    % initial guess of initial value eta
    eta = [kc1_jig;0;0];

    % initial cost function value
%     [wv,pmat,Rmat] = fn_intgEP_v1_2layers_kconst(eta,kc1_jig,kc2_jig,z_crit,theta0,ds,N,B,Binv);
    scalef0 = 1;
    Cval = costfn_opt_v2_singlebend_2layers_kconst(eta,kc1_jig,kc2_jig,...
        data_jig,s_m,0,ds,theta0,z_crit,N,B,Binv,scalef0);
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
        'MaxFunEvals',10000,'Display','notify','MaxIter',maxiter,'MaxFunEvals',100000);

    %psoptions = psoptimset(psoldopts,'Display','iter','MaxIter',maxiter);
    % gaset = saoptimset;
    % gaset = saoptimset(gaset,'Display','iter','TimeLimit',timelimit,'TolFun',tolfun);

    disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

    x0 = eta; % initial value
    LB = [-0.01*ones(3,1)]; % lower bound
    UB = [0.01*ones(3,1)]; % upper bound

    [x,fval,exitflag] = fmincon(@(x) costfn_opt_v2_singlebend_2layers_kconst(x,...
        kc1_jig,kc2_jig,data_jig,s_m,0,ds,theta0,z_crit,N,B,Binv,scalef),...
        x0,[],[],[],[],LB,UB,[],options);

    t_f = toc;

    % x
    % fval
    % exitflag

    disp(['%>>>> optimization ends: time = ',num2str(t_f),' [sec]'])

    % final configuration
    disp(['%>> w_init = '])
    disp(num2str(x'))

    [wv_f,pmat,Rmat] = fn_intgEP_v1_2layers(x,kc1_jig,kc2_jig,z_crit,theta0,0,ds,N,B,Binv);
    
    % save
    wmat_tot_jig(:,:,ii) = wv_f;
    
    %% beam data
    disp(['%%%%% beam data %%%%%'])

    data_mat = Data_raw_beam(:,ii);
    w_m1 = [data_mat(1);data_mat(2);0];
    w_m2 = [data_mat(3);data_mat(4);0];
    w_m3 = [data_mat(5);data_mat(6);0];
    data_beam = [w_m1,w_m2,w_m3]*1e-3; % in 1/mm

    data_beam_tot(:,:,ii) = data_beam;
    
    % initial guess of initial value eta
    eta = [kc1_beam;0;0];
    
    % initial cost function value
%     [wv,pmat,Rmat] = fn_intgEP_v1_2layers(eta,kc1_beam,kc2_beam,z_crit,theta0,0,ds,N,B,Binv);
    scalef0 = 1;
    Cval = costfn_opt_v2_singlebend_2layers_kconst(eta,kc1_beam,kc2_beam,...
        data_beam,s_m,0,ds,theta0,z_crit,N,B,Binv,scalef0);
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
        'MaxFunEvals',10000,'Display','notify','MaxIter',maxiter,'MaxFunEvals',100000);

    %psoptions = psoptimset(psoldopts,'Display','iter','MaxIter',maxiter);
    % gaset = saoptimset;
    % gaset = saoptimset(gaset,'Display','iter','TimeLimit',timelimit,'TolFun',tolfun);

    disp(['%>>>> Optimization (fmincon) is being executed... Please wait...'])

    x0 = eta; % initial value
    LB = [-0.01*ones(3,1)]; % lower bound
    UB = [0.01*ones(3,1)]; % upper bound

    [x,fval,exitflag] = fmincon(@(x) costfn_opt_v2_singlebend_2layers_kconst(x,...
        kc1_beam,kc2_beam,data_beam,s_m,0,ds,theta0,z_crit,N,B,Binv,scalef),...
        x0,[],[],[],[],LB,UB,[],options);

    t_f = toc;

    % x
    % fval
    % exitflag

    disp(['%>>>> optimization ends: time = ',num2str(t_f),' [sec]'])

    % final configuration
    disp(['%>> w_init = '])
    disp(num2str(x'))

    [wv_f,pmat,Rmat] = fn_intgEP_v1_2layers(x,kc1_beam,kc2_beam,z_crit,theta0,0,ds,N,B,Binv);

    % save
    wmat_tot_beam(:,:,ii) = wv_f;

    disp([' '])
end % end of total loop

%% extract needle trajectories from image analysis data for comparison
% % note on file names
% pos_img = xlsread('T2009.xlsx','06','D:E'); % , n180 case, number 6
% pos_img = xlsread('T0919.xlsx','T0919_06','D:E'); % number 6, sign should change
% pos_img_tot = zeros(83,2,tot_num_exp);
pos_img_tot = cell(tot_num_exp,1);
for jj = 1:tot_num_exp
    if jj <= 9
        sheet_num = ['0',num2str(jj)];
    elseif jj == 10
        sheet_num = ['10'];
    end
%     pos_img = xlsread('T0919.xlsx',sheet_num,'C:D');
    pos_img = readmatrix('T0919.xlsx','Sheet',sheet_num,'Range','C:D');
    pos_img_tot{jj} = pos_img;
end
% last point is out-of-range

%% optimization for initial rotation angle
% initialize
theta0_jig = zeros(1,tot_num_exp);
theta0_beam = zeros(1,tot_num_exp);

pmat_tot_jig = zeros(3,N,tot_num_exp);
pmat_tot_beam = zeros(3,N,tot_num_exp);

% initial theta0
theta0 = 1.0*180/pi; % 1 deg

for ii = 1:tot_num_exp
    pos_img = pos_img_tot{ii};

    %% jig calibration results
    wvec = wmat_tot_jig(:,:,ii);
    
    % initial cost function value
    cost_th0 = costfn_theta0(theta0,wvec,pos_img,ds,1);
    
    scalef = 1/cost_th0;

    % optimization
    tic;
    maxiter = 1000;
    timelimit = 100;
    tolfun = 1e-6;%1e-4;

    oldopts = optimset('fmincon');
    Tol = 1e-14*10^(ceil(log10(scalef)));

    options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
        'MaxFunEvals',10000,'Display','notify');

    disp(['%>>>> Optimization (fmincon) on theta0 is being executed... Please wait...'])

    th0 = theta0; % initial value
    LB = 0; % lower bound
    UB = pi/10; % upper bound

    [th,fval,exitflag] = fmincon(@(th) costfn_theta0(th,wvec,pos_img,ds,scalef),...
        th0,[],[],[],[],LB,UB,[],options);
    
    tf = toc;
    disp(['%>>>> optimization ends: time = ',num2str(tf),' [sec]'])
    
    % save
    theta0_jig(ii) = th;
    
    % configuration
    Rmat = zeros(3,3,N);
    Rmat(:,:,1) = Rot_x(th);%eye(3);
    pmat = zeros(3,N);
    for i = 2:N
        % orientation
        W = matr(1/2*(wvec(:,i-1) + wvec(:,i)));
        Rmat(:,:,i) = Rmat(:,:,i-1)*expm(ds*W);

        % position
        e3vec = squeeze(Rmat(:,3,1:i));
        if i == 2
            pmat(:,i) = pmat(:,i-1) + squeeze(Rmat(:,3,i))*ds;
        else
            pmat(:,i) = Simpson_vec_int(e3vec,ds);
        end        
    end

    pmat_tot_jig(:,:,ii) = pmat;

    %% beam calibration results
    wvec = wmat_tot_beam(:,:,ii);
    
    % initial cost function value
    cost_th0 = costfn_theta0(theta0,wvec,pos_img,ds,1);
    
    scalef = 1/cost_th0;

    % optimization
    tic;
    maxiter = 1000;
    timelimit = 100;
    tolfun = 1e-6;%1e-4;

    oldopts = optimset('fmincon');
    Tol = 1e-14*10^(ceil(log10(scalef)));

    options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
        'MaxFunEvals',10000,'Display','notify');

    disp(['%>>>> Optimization (fmincon) on theta0 is being executed... Please wait...'])

    th0 = theta0; % initial value
    LB = 0; % lower bound
    UB = pi/10; % upper bound

    [th,fval,exitflag] = fmincon(@(th) costfn_theta0(th,wvec,pos_img,ds,scalef),...
        th0,[],[],[],[],LB,UB,[],options);
    
    tf = toc;
    disp(['%>>>> optimization ends: time = ',num2str(tf),' [sec]'])
    
    % save
    theta0_beam(ii) = th;
    
    % configuration
    Rmat = zeros(3,3,N);
    Rmat(:,:,1) = Rot_x(th);%eye(3);
    pmat = zeros(3,N);
    for i = 2:N
        % orientation
        W = matr(1/2*(wvec(:,i-1) + wvec(:,i)));
        Rmat(:,:,i) = Rmat(:,:,i-1)*expm(ds*W);

        % position
        e3vec = squeeze(Rmat(:,3,1:i));
        if i == 2
            pmat(:,i) = pmat(:,i-1) + squeeze(Rmat(:,3,i))*ds;
        else
            pmat(:,i) = Simpson_vec_int(e3vec,ds);
        end        
    end

    pmat_tot_beam(:,:,ii) = pmat;
end

%% theta0 distribution
disp(' ')
disp(['average \theta_0 (jig) = ',num2str(mean(theta0_jig)*180/pi),' [deg]']);
disp(['std \theta_0 (jig) = ',num2str(std(theta0_jig)*180/pi),' [deg]']);
disp(' ')
disp(['average \theta_0 (beam) = ',num2str(mean(theta0_beam)*180/pi),' [deg]']);
disp(['std \theta_0 (beam) = ',num2str(std(theta0_beam)*180/pi),' [deg]']);
disp(' ')

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
    
    tip_defl1 = abs(pos_img_t_jig(end,2) + pmat_tot_jig(2,end,jj));
    tip_defl2 = abs(pos_img_t_beam(end,2) + pmat_tot_beam(2,end,jj));
    
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
    
    y_data_jig = -pos_img_t_jig(:,2);
    z_data_jig = pos_img_t_jig(:,1);
    
    y_data_beam = -pos_img_t_beam(:,2);
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
    plot(pos_img_true_jig{i_plot}(1:2:end,1),-pos_img_true_jig{i_plot}(1:2:end,2),'ko','MarkerSize',8)
    plot(pos_img_true_beam{i_plot}(1:2:end,1),-pos_img_true_beam{i_plot}(1:2:end,2),'ko','MarkerSize',8)
    hold off
    grid on
    axis equal
    axis([0 100 -12 10])
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

% %pos_img = xlsread('T2009.xlsx','06','D:E'); % , n180 case, number 6
% pos_img = xlsread('T0919.xlsx','T0919_06','D:E'); % number 6, sign should change
% jj = 6;
% 
% tip_defl1 = abs(-pos_img(end,2) - pmat_tot_jig(2,end,jj));
% tip_defl2 = abs(-pos_img(end,2) - pmat_tot_beam(2,end,jj));
% 
% disp(['%%%%%%%%%%%% number of experiment = ',num2str(jj),' %%%%%%%%%%%'])
% disp(['difference of tip deflection = ',num2str(tip_defl1),' [mm] (jig)'])
% disp(['difference of tip deflection = ',num2str(tip_defl2),' [mm] (beam)'])
% disp([' '])
% 
% y_data = -pos_img(:,2);
% z_data = pos_img(:,1);
% 
% pos_jig = pmat_tot_jig(:,:,jj);
% pos_beam = pmat_tot_beam(:,:,jj);
% 
% RMSE_jig = extract_data_pos(y_data,z_data,pos_jig);
% RMSE_beam = extract_data_pos(y_data,z_data,pos_beam);
% 
% disp(['%%%%%%%%%%%% number of experiment = ',num2str(jj),' %%%%%%%%%%%'])
% disp(['RMSE = ',num2str(RMSE_jig),' [mm] (jig)'])
% disp(['RMSE = ',num2str(RMSE_beam),' [mm] (beam)'])
% disp([' '])

% i_plot = 6;
% figure;
% subplot(3,1,1)
% plot(s,wmat_tot_jig(1,:,i_plot),'k-','LineWidth',2)
% hold on
% plot(s,wmat_tot_beam(1,:,i_plot),'b-','LineWidth',2)
% plot(s(s_m),data_jig_tot(1,:,i_plot),'ko','MarkerSize',10,'LineWidth',2)
% plot(s(s_m),data_beam_tot(1,:,i_plot),'bs','MarkerSize',10,'LineWidth',2)
% hold off
% grid on
% ylabel('{\omega}_x [1/mm]')
% xlabel('Arclength [mm]')        
% legend('model (jig)','model (beam)','measured (jig)','measured (beam)')
% 
% subplot(3,1,2)
% plot(s,wmat_tot_jig(2,:,i_plot),'k-','LineWidth',2)
% hold on
% plot(s,wmat_tot_beam(2,:,i_plot),'b-','LineWidth',2)
% plot(s(s_m),data_jig_tot(2,:,i_plot),'ko','MarkerSize',10,'LineWidth',2)
% plot(s(s_m),data_beam_tot(2,:,i_plot),'bs','MarkerSize',10,'LineWidth',2)
% hold off
% grid on
% ylabel('{\omega}_y [1/mm]')
% xlabel('Arclength [mm]')
% legend('model (jig)','model (beam)','measured (jig)','measured (beam)')
% 
% subplot(3,1,3)
% plot(squeeze(pmat_tot_jig(3,:,i_plot)),squeeze(pmat_tot_jig(2,:,i_plot)),'k-','LineWidth',2)
% hold on
% plot(squeeze(pmat_tot_beam(3,:,i_plot)),squeeze(pmat_tot_beam(2,:,i_plot)),'b-','LineWidth',2)
% plot(pos_img(1:2:end,1),-pos_img(1:2:end,2),'ko','MarkerSize',8)
% hold off
% grid on
% axis equal
% axis([0 100 -12 8])
% ylabel('Y [mm]')
% xlabel('Z [mm]')
% legend('model (jig)','model (beam)','image')




