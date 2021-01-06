%% fn_Workspace_2layer_data_Dimitri.m
%
%
% used for the simulation of workspace for data-based double layer needle insertion.
%
% - written by Jin Seob Kim
% - edited  by Dimitri Lezcano

function pmat_total = fn_Workspace_2layer_data_Dimitri(exp_num, cal_type, save_bool)
    global p L lengths ds theta0 B Binv w_init_ref s_crit N_arclengths z_crit 

    % saving set-up
    directory = "Dimitri/Data/";
    directory = directory + "Archive/2-Layer/TipAError_Cost_Data/";
    file_datatable = directory + "Double_Layer_Data.mat"; % file for expmt. data table

    directory = directory + "kc_p_w_init_q/";
    directory = directory + sprintf("exp_%02d/", exp_num);
    file_save = directory + "kc_doublelayer_workspace_data";

    file_save = file_save + "_" + cal_type;
    disp(file_save);

    %% load the data table and the associated data
    exp_data_tbl = load(file_datatable, 'l2_data').l2_data;
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
    ds = 0.5; % in mm, 1 mm seems fine as well
    L = 90; % in mm
    s = [0:ds:L];
    N = length(s);

    %% intrinsic curvature
    % first reference value when L = 90 mm
    if strcmp(data_row.type, "jig")
        kc_s = 0.0017;
        kc_h = 0.0027;
        theta0 = 1.1287e-09 * pi/180;

    elseif strcmp(data_row.type, "beam")
        kc_s = 0.0020;
        kc_h = 0.0032;
        theta0 = 0.0027919 * pi/180;
    end


    % data
    w_init_ref = data_row.w_init_ref' ; % need to change

    %% generate arclength array
    lengths = 90:15:150;
    N_arclengths = floor(lengths/ds) + 1;

    %% determine s_crit (from 90 mm homogeneous single bend insertion)
    z_crit = 45;
    [~,pos1] = fn_intgEP_v1_1layer([kc_s;0;0],kc_s,0,0,ds,N,B,Binv);
    ix_crit_v = find(abs(pos1(3,:) - z_crit) <= ds/2); 
    ix_crit = ix_crit_v(1);
    s_crit = s(ix_crit);

    %% shape set-up as an array
    p = 0.5860;

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

    cost_function = @(q) cost_fn(q, kc_s, kc_h);

    q0 = p;
    [q, min_err, exitflag, output] = fmincon(cost_function, q0, [],[],[],[],low_bound,up_bound,[],options);

    %% Generate rotated shapes
    L_pred = 150; N_pred = L_pred/ds + 1;
    pmat_total = cell(1, length(theta_rot) + 1);

    % generate the reference L_ref insertion
    [~, pmat_total{1}, ~] = fn_intgEP_v2_2layers(w_init_ref,kc_s,kc_h,z_crit,theta0,0,ds,N_arclengths(1),B,Binv);
                            
    % generate the rotated predictions
    for i = 1:length(theta_rot)
        th_rot = theta_rot(i);

        % generate the rotated shape
        pmat_total{i+1} = integrate_twist(kc_s, kc_h, L_pred, N_pred, q, th_rot);

    end

    %% pair the circles for conical look
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

    % patch for boundary of tissue
    [~, idx] = min(abs(pmat_i(3,:)-z_crit)); % find the z_crit index
    S1 = [z_crit,        z_crit;
          pmat_i(1,idx) - 5, pmat_i(1,idx) + 5;
          pmat_i(1,idx) - 5 , pmat_i(2,idx) + 1];
    S2 = S1;
    S2(2,:) = S1(2,end:-1:1);
    S = [S1(:,1) S2(:,1) S1(:,2) S2(:,2)];

    patch(S(1,:), S(2,:), S(3,:), 'r', 'DisplayName', 'Tissue Boundary');
    
    % plotting options
    set(gcf,'units', 'normalized', 'position', [1/10, 1/10, 6/8, 6/8]);
    xlabel('z [mm]'); ylabel('x [mm]'); zlabel('y [mm]')
    grid on; axis equal; 
    view([53.966666, 8.891815]);
    hold off;

    %% saving
   
    if save_bool
        % figure
        savefig(f1, file_save + '_workspace_3d.fig');
        fprintf('Saved figure: %s\n', file_save + '_workspace_3d.fig');
        saveas(f1, file_save + '_workspace_3d.png');
        fprintf('Saved image: %s\n', file_save + '_workspace_3d.png');

    end
    

end
%% Functions
% cost function for optimization
function cost = cost_fn(q, kc_s, kc_h)    
    global p L lengths ds theta0 N_arclengths B Binv w_init_ref  z_crit
    
    pmat_total = cell(1,length(lengths));
    
    % needle trajectories for all lengths
    for i = 1:length(lengths)
        Li = lengths(i);
        Ni = N_arclengths(i);
        
        % generate the shape
        kc_s_i = kc_s * (L/Li)^p;
        kc_h_i = kc_h * (L/Li)^p;
        w_init = w_init_ref * (L/Li)^q;
        [~,pmati,~] = fn_intgEP_v3_2layers(w_init,kc_s_i,kc_h_i,z_crit,theta0,0,ds,Ni,B,Binv);
        
        % add the shape to a cell
        pmat_total{i} = pmati;
        
    end
    
    cost = TipAerror_Dimitri(pmat_total, lengths);
    
end

% function for integrating using a rotation at reference L
function pmat = integrate_twist(kc1, kc2, Li, Ni, q, theta_rot)
    global w_init_ref N_arclengths L p ds B Binv theta0 s_crit
    
    % parameter scaling
    kc1_i = kc1*(L/Li)^p; % scale the kappa_c value
    kc2_i = kc2*(L/Li)^p; % scale the kappa_c value
    w_init = w_init_ref * (L/Li)^q;
    
    % arclength
    s = 0:ds:Li;
    [~,ix_crit] = min(abs(s - s_crit));
    s1 = s(1:ix_crit);
    s2 = s(ix_crit+1:Ni);
    
    % intrinsic curvature
    k0 = zeros(1,Ni);
    k0(1:ix_crit) = kc1_i*(s_crit - s1).^2/Li^2 + kc2_i*(1 - s_crit/Li)*(1 + s_crit/Li - 2*s1/Li);
    k0(ix_crit+1:Ni) = kc2_i*(1 - s2/Li).^2;
    w0 = [k0;zeros(1,Ni);zeros(1,Ni)];

    k0prime = zeros(1,Ni);
    k0prime(1:ix_crit) = -2*kc1_i/Li^2*(s_crit - s1) - 2*kc2_i/Li*(1 - s_crit/Li);
    k0prime(ix_crit+1:Ni) = -2*kc2_i/Li*(1 - s2/Li);
    w0prime = [k0prime;zeros(1,Ni);zeros(1,Ni)];

    % Rotate using rotation z matriz 
    Rz = Rot_z(theta_rot);
    w0(:,N_arclengths(1)+1:end) = Rz*w0(:,N_arclengths(1)+1:end);
    w0prime(:,N_arclengths(1)+1:end) = Rz*w0prime(:,N_arclengths(1)+1:end);

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