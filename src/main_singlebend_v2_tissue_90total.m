%% main_singlebend_v2_tissue_90total.m
%
% old comment: use needleshape_fmincon_v2.m and data 09-02-2016
% this is the one used for ICRA 2017.
% single-bend, single homogeneous layer case
%
% - written by Jin Seob (Jesse) Kim

clear all

set(0,'DefaultAxesFontSize',20);

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


%% load data
% New data on 09-06-2016 hard phantom [1/m] (90 mm insertion, in tissue)
tot_num_exp = 10;% # of experiments 

% Begin the results data table
results_tbl = table('Size', [2*tot_num_exp, 8], 'VariableTypes', ["int8", "string", "double", "double", "double", "double", "double", "double"], ...
    'VariableNames', ["exp_num", "cal_type", "kc", "w_init1", "w_init2", "w_init3", "tip_err", "rmse"]);
results_tbl = mergevars(results_tbl, 4:6, 'NewVariableName', 'w_init_ref');

% load the data
Data_raw_tot = zeros(6,20,tot_num_exp);
for i_exp = 1:tot_num_exp %SLOW HERE, OLD FUNCTION USED - DIMITRI % FIXED - DIMITRI
    disp("Processing Experiment #"); disp(i_exp)
    sheet_num = ['T',num2str(i_exp)];
%     data_tmp = xlsread('Results_insertion 0906_HardPhantom.xlsx',sheet_num,'C:V');
%     Data_raw_tot(:,:,i_exp) = data_tmp(3:8,:);
    % Faster processing time
%     data_tmp = readmatrix('real_data_hardphantom.xlsx','Sheet',sheet_num,'Range','C:V');
    data_tmp = readmatrix('New_Frame_Results_Insertion 0902.xlsx','Sheet',sheet_num,'Range','C:V');
    
    Data_raw_tot(:,:,i_exp) = data_tmp(1:6,:);
end

% % New data on 09-02-2016 [1/m] (90 mm insertion, in tissue)
% tot_num_exp = 10;
% 
% Data_raw_tot = zeros(6,20,tot_num_exp);
% for i_exp = 1:tot_num_exp
%     sheet_num = ['T',num2str(i_exp)];
%     data_tmp = xlsread('New_Frame_Results_insertion 0902.xlsx',sheet_num,'C:V');
%     Data_raw_tot(:,:,i_exp) = data_tmp(3:8,:);
% end

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
theta0 = 1.0;%0.7;, where to choose this? 
theta0 = theta0*pi/180;

%% total loop
kappa_c_mat = zeros(tot_num_exp,2);

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
    
    % parameters on intrinsic curvature (initial values)
    kc = data_jig(1,1)*(1+0.7);

    % initial configuration of a needle
    % initial guess of initial value eta
    eta = zeros(4,1);
    eta(1:3) = [kc;0;0];
    eta(4) = kc;

    % intrinsic curvature
    k0 = kc*(1 - s/L).^2;
    w0 = [k0;zeros(1,N);zeros(1,N)];

    k0prime = -2*kc/L*(1 - s/L);
    w0prime = [k0prime;zeros(1,N);zeros(1,N)];

    % initial cost value
    w_init = eta(1:3);
    wv = fn_intgEP_v1(w_init,w0,w0prime,0,ds,N,B,Binv);
    scalef0 = 1;
    Cval = costfn_opt_v2_singlebend(eta,data_jig,s_m,ds,N,B,Binv,scalef0);
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
    LB = [-0.01*ones(3,1);0]; % lower bound
    UB = [0.01*ones(3,1);0.01]; % upper bound

    [x,fval,exitflag] = fmincon(@(x) costfn_opt_v2_singlebend(x,data_jig,s_m,ds,N,B,Binv,scalef),...
        x0,[],[],[],[],LB,UB,[],options);

    t_f = toc;

%     x
%     fval
%     exitflag

    disp(['%>>>> optimization ends: time = ',num2str(t_f),' [sec]'])

    % final configuration
    kc_f = x(4);
    kappa_c_mat(ii,1) = kc_f;
    
    % append to the table
    results_tbl(2*ii-1,1:4) = {ii, "jig", kc_f, reshape(x(1:3), 1, 3)};
    
    disp(['%>>>> kappa_c = ',num2str(kc_f)])

    disp(['%>> w_init = '])
    disp(num2str(x(1:3)))

    k0_f = kc_f*(1 - s/L).^2;
    w0_f = [k0_f;zeros(1,N);zeros(1,N)];

    k0prime_f = -2*kc_f/L*(1 - s/L);
    w0prime_f = [k0prime_f;zeros(1,N);zeros(1,N)];

    wv_f = fn_intgEP_v1(x(1:3),w0_f,w0prime_f,0,ds,N,B,Binv);

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
    w_m1 = [data_mat(1);data_mat(2);0];
    w_m2 = [data_mat(3);data_mat(4);0];
    w_m3 = [data_mat(5);data_mat(6);0];
    data_beam = [w_m1,w_m2,w_m3]*1e-3; % in 1/mm

    data_beam_tot(:,:,ii) = data_beam;

    % parameters on intrinsic curvature (initial values)
    kc = data_beam(1,1)*(1+0.7);

    % initial configuration of a needle
    % initial guess of initial value eta
    eta = zeros(4,1);
    eta(1:3) = [kc;0;0];
    eta(4) = kc;

    % intrinsic curvature
    k0 = kc*(1 - s/L).^2;
    w0 = [k0;zeros(1,N);zeros(1,N)];

    k0prime = -2*kc/L*(1 - s/L);
    w0prime = [k0prime;zeros(1,N);zeros(1,N)];

    % initial cost value
    w_init = eta(1:3);
    wv = fn_intgEP_v1(w_init,w0,w0prime,0,ds,N,B,Binv);
    scalef0 = 1;
    Cval = costfn_opt_v2_singlebend(eta,data_beam,s_m,ds,N,B,Binv,scalef0);
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
    LB = [-0.01*ones(3,1);0]; % lower bound
    UB = [0.01*ones(3,1);0.01]; % upper bound

    [x,fval,exitflag] = fmincon(@(x) costfn_opt_v2_singlebend(x,data_beam,s_m,ds,N,B,Binv,scalef),...
        x0,[],[],[],[],LB,UB,[],options);

    t_f = toc;

    % x
    % fval
    % exitflag

    disp(['%>>>> optimization ends: time = ',num2str(t_f),' [sec]'])

    % final configuration
    kc_f = x(4);
    kappa_c_mat(ii,2) = kc_f;
    
    disp(['%>>>> kappa_c = ',num2str(kc_f)])
     
    % append to the table
    results_tbl(2*ii,1:4) = {ii, "beam", kc_f, reshape(x(1:3), 1, 3)};

    disp(['%>> w_init = '])
    disp(num2str(x(1:3)))

    k0_f = kc_f*(1 - s/L).^2;
    w0_f = [k0_f;zeros(1,N);zeros(1,N)];

    k0prime_f = -2*kc_f/L*(1 - s/L);
    w0prime_f = [k0prime_f;zeros(1,N);zeros(1,N)];

    wv_f = fn_intgEP_v1(x(1:3),w0_f,w0prime_f,0,ds,N,B,Binv);

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
% % soft phantom
% pos_img_tot = zeros(83,2,tot_num_exp);
% for jj = 1:tot_num_exp
%     if jj <= 9
%         sheet_num = ['0902_0',num2str(jj)];
%     elseif jj == 10
%         sheet_num = ['0902_10'];
%     end
%     pos_img = xlsread('Results_T0902.xlsx',sheet_num,'D:E');
%     pos_img_tot(:,:,jj) = pos_img;
% end

% hard phantom
pos_img_tot = zeros(83,2,tot_num_exp);
for jj = 1:tot_num_exp
    if jj <= 9
        sheet_num = ['0',num2str(jj)];
    elseif jj == 10
        sheet_num = ['10'];
    end
%     pos_img = xlsread('T0906_HardPhantom.xlsx',sheet_num,'C:D');
%     pos_img = readmatrix('T0906_HardPhantom.xlsx','Sheet',sheet_num,'Range','C:D');
    pos_img = readmatrix('Results_T0902.xlsx','Sheet',['0902_',sheet_num],'Range','C:D');
    num_rows_img = length(pos_img);
    pos_img_tot(1:num_rows_img,:,jj) = pos_img;
end

%% tip deflection calculation
tip_defl_jig = zeros(tot_num_exp,1);
tip_defl_beam = zeros(tot_num_exp,1);

for jj = 1:tot_num_exp
    tip_defl1 = abs(pos_img_tot(end,2,jj) + pmat_tot_jig(2,end,jj));
    tip_defl2 = abs(pos_img_tot(end,2,jj) + pmat_tot_beam(2,end,jj));
    
    tip_defl_jig(jj) = tip_defl1;
    tip_defl_beam(jj) = tip_defl2;
    
    % append to table
    results_tbl{2*jj-1, 5} = tip_defl_jig(jj);
    results_tbl{2*jj, 5} = tip_defl_beam(jj);

    disp(['%%%%%%%%%%%% number of experiment = ',num2str(jj),' %%%%%%%%%%%'])
    disp(['difference of tip deflection = ',num2str(tip_defl1),' [mm] (jig)'])
    disp(['difference of tip deflection = ',num2str(tip_defl2),' [mm] (beam)'])
    disp([' '])
end

%% root mean square error of the trajectories
RMSE_jig = zeros(tot_num_exp,1);
RMSE_beam = zeros(tot_num_exp,1);

for jj = 1:tot_num_exp
    y_data = -pos_img_tot(:,2,jj);
    z_data = pos_img_tot(:,1,jj);
    
    pos_jig = pmat_tot_jig(:,:,jj);
    pos_beam = pmat_tot_beam(:,:,jj);
    
    RMSE_jig(jj) = extract_data_pos(y_data(1:end-1),z_data(1:end-1),pos_jig);
    RMSE_beam(jj) = extract_data_pos(y_data(1:end-1),z_data(1:end-1),pos_beam);
    
    % append to table
    results_tbl{2*jj-1, 6} = RMSE_jig(jj);
    results_tbl{2*jj, 6} = RMSE_beam(jj);

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
    plot(pos_img_tot(1:2:end,1,i_plot),-pos_img_tot(1:2:end,2,i_plot),'ko','MarkerSize',8)
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


