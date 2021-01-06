%% fn_Workspace_1layer_data_Dimitri.m
%
%
% used for the simulation of workspace for data-based single bend needle insertion.
%
% - written by Jin Seob Kim
% - edited  by Dimitri Lezcano

function pmat_total = fn_Workspace_1layer_data_Dimitri(exp_num, cal_type, save_bool, exp_type)
    global p L lengths ds theta0 N_arclengths B Binv w_init_ref 

    % saving set-up
    directory = "Dimitri/Data/";
    directory = directory + "Archive/1-Layer/TipAError_Cost_Data/";
    if strcmpi(exp_type, "soft")
        file_datatable = directory + "Single_Layer_Data_Soft.mat";
    else
        file_datatable = directory + "Single_Layer_Data.mat"; % file for expmt. data table
        exp_type = "";
    end
    
    directory = directory + "kc_p_w_init_q/";
    directory = directory + sprintf("exp_%02d/", exp_num);
    file_save = directory + "kc_singlelayer_workspace_data";
    
    if strcmpi(exp_type, "soft")
        file_save = file_save + "_" + lower(exp_type);
    end
  
    file_save = file_save + "_" + cal_type;
    disp(file_save);

    %% load the data table and the associated data
    exp_data_tbl = load(file_datatable, 'l1_data').l1_data;
    data_row = exp_data_tbl(strcmp( ...
        exp_data_tbl.type, cal_type) & exp_data_tbl.exp_num == exp_num,:);

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


    %% intrinsic curvature
    % first reference value when L = 90 mm
    kc_s = data_row.kc;
    p = 0.5860;

    % data
    w_init_ref = reshape(data_row.w_init_ref,3,1);
    theta0 = 1.0*pi/180;

    %% generate arclength array
    ds = 0.5; % in mm
    L1 = 90; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
    L2 = 105; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
    L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
    L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
    L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

    lengths = [L1, L2, L3, L4, L5];
    N_arclengths = lengths/ds + 1;

    %% rotation values as an array
    theta_rot = linspace(0, 2*pi, 50);

    %=========================================================================
    %% Determine optimal q
    low_bound = 0.0;
    up_bound = 1.0;

    Tol = 1e-14;
    opts = optimset('fmincon');
    options = optimset(opts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
            'MaxFunEvals',10000,'Display','notify');

    cost_function = @(q) cost_fn(q, kc_s);

    q0 = p;
    [q, min_err, exitflag, output] = fmincon(cost_function, q0, [],[],[],[],low_bound,up_bound,[],options);

    %% Generate rotated shapes
    L_pred = 150; N_pred = L_pred/ds + 1;
    pmat_total = cell(1, length(theta_rot) + 1);

    % generate the reference L_ref insertion
    [~, pmat_total{1}, ~] = fn_intgEP_v1_1layer(w_init_ref,kc_s,theta0,0,ds,N_arclengths(1),B,Binv);

    % generate the rotated predictions
    for i = 1:length(theta_rot)
        th_rot = theta_rot(i);

        % generate the rotated shape
        pmat_total{i+1} = integrate_twist(kc_s, L_pred, N_pred, q, th_rot);

    end

    %% find the circles for cone fitting
    s_pred = L:ds:L_pred;
    circle_fit = cell(1, length(s_pred));

    % find the circle fittings
    for i = 1:length(circle_fit)
        points = zeros(3,length(pmat_total)-1);
        for j = 2:length(pmat_total)
            pmat_j = pmat_total{j};
            points(:,j-1) = pmat_j(:,N_arclengths(1) + i - 1);

        end
        circle_fit{i} = points;
    end

    %% position plots
    f1 = figure(1);
    for i = length(circle_fit):-1:2
        % get the position
    %     pmat_i = pmat_total{i};
        circle_i = circle_fit{i};

        plot3(circle_i(3,:), circle_i(1,:), circle_i(2,:), 'k-o', 'linewidth', 2,'markersize',.1);
        hold on;

    end

    % plot the 90 mm case
    pmat_i = pmat_total{1};
    plot3(pmat_i(3,:), pmat_i(1,:), pmat_i(2,:), 'k', 'linewidth', 2);

    % % plot the first scenario th_rot = 0
    % pmat_i = pmat_total{2};
    % plot3(pmat_i(3,:), pmat_i(1,:), pmat_i(2,:), 'k--', 'linewidth', 1);

    xlabel('z [mm]'); ylabel('x [mm]'); zlabel('y [mm]')
    grid on; axis equal; 
    view([53.966666, 8.891815]);
    hold off;

    %% saving
    msg = "";
    msg = msg + sprintf('kc(90 mm), %.6f\n', kc_s);
    msg = msg + sprintf('w_init(90 mm), [ %.5f; %.5f; %.5f ]\n', w_init_ref);
    msg = msg + sprintf('p, %.5f\n', p);
    msg = msg + sprintf('best q, %.5f\n', q);
    msg = msg + sprintf('TipAError (mm^2), %.3f\n', min_err);
    
    if save_bool
        % figure
        savefig(f1, file_save + '_workspace_3d.fig');
        fprintf('Saved figure: %s\n', file_save + '_workspace_3d.fig');
        saveas(f1, file_save + '_workspace_3d.png');
        fprintf('Saved image: %s\n', file_save + '_workspace_3d.png');

        % csv file
        fid = fopen(file_save + '_err.csv', 'w');
        fprintf(fid, msg);
        fprintf('Wrote file: %s\n', file_save + '_err.csv'); 
        fclose(fid);
    end
    
    disp(msg);
    
end
%% Functions
% cost function for optimization
function cost = cost_fn(q, kc)    
    global p L lengths ds theta0 N_arclengths B Binv w_init_ref 
    
    pmat_total = cell(1,length(lengths));
    
    % needle trajectories for all lengths
    for i = 1:length(lengths)
        Li = lengths(i);
        Ni = N_arclengths(i);
        
        % generate the shape
        kc_i = kc * (L/Li)^p;
        w_init = w_init_ref * (L/Li)^q;
        [~,pmati,~] = fn_intgEP_v1_1layer(w_init,kc_i,theta0,0,ds,Ni,B,Binv);

        % add the shape to a cell
        pmat_total{i} = pmati;
        
    end
    
    cost = TipAerror_Dimitri(pmat_total, lengths);
    
end

% function for integrating using a rotation at reference L
function pmat = integrate_twist(kc, Li, Ni, q, theta_rot)
    global w_init_ref N_arclengths L p ds B Binv theta0
    
    % parameter scaling
    kc_i = kc*(L/Li)^p; % scale the kappa_c value
    w_init = w_init_ref * (L/Li)^q;
    
    % arclength
    s = 0:ds:Li;
    
    % intrinsic curvature
    k0 = kc_i*(1 - s/Li).^2;
    w0 = [k0;zeros(1,Ni);zeros(1,Ni)];

    k0prime = -2*kc_i/Li*(1 - s/Li);
    w0prime = [k0prime;zeros(1,Ni);zeros(1,Ni)];

    % Rotate using rotation z matriz 
    Rz = Rot_z(theta_rot);
    w0(:,N_arclengths(1)+1:end) = Rz*w0(:,N_arclengths(1)+1:end);
    w0prime(:,N_arclengths(1)+1:end) = Rz*w0prime(:,N_arclengths(1)+1:end);

    % integrate EP equation
    [wv,~,~] = fn_intgEP_w0_Dimitri(w_init, w0, w0prime,theta0,0,ds,Ni,B,Binv);
    
    % rotate the omega vector
    wv_rot = wv;
    wv_rot(:,N_arclengths(1)+1:end) = Rz * wv_rot(:,N_arclengths(1)+1:end);
    
    % integrate the omega vector to get positions
    [pmat, ~] = intg_MatDiffEq(wv_rot);
end

% integrates wv here essentially
function [pmat, Rmat] = intg_MatDiffEq(wv)

global ds theta0

N = size(wv,2);
Rmat = zeros(3,3,N);
Rmat(:,:,1) = Rot_x(theta0);

pmat = zeros(3,N);
for i = 2:N
    % orientation
    W = matr(1/2*(wv(:,i-1) + wv(:,i)));
    Rmat(:,:,i) = Rmat(:,:,i-1)*expm(ds*W);

    % position
    e3vec = squeeze(Rmat(:,3,1:i));
    if i == 2
        pmat(:,i) = pmat(:,i-1) + squeeze(Rmat(:,3,i))*ds;
    else
        pmat(:,i) = Simpson_vec_int(e3vec,ds);
    end        
end

end