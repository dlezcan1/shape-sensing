%% fn_ErrorComp_q_2layer_data_Dimitri.m
%
% from ideal_singlebendinsertion_2020_0413.m
%
% used for the simulation of data-based single bend needle insertion.
%
% - written by Dimitri Lezcano

function [pmat_total, q_best, min_err] = fn_ErrorComp_q_2layer_data_Dimitri(exp_num, cal_type, save_bool)
    global p L lengths ds theta0 B Binv w_init_ref s_crit N_arclengths z_crit

    % saving set-up
    directory = "Dimitri/Data/";
    directory = directory + "Archive/2-Layer/TipAError_Cost_Data/";
    file_datatable = directory + "Double_Layer_Data.mat"; % file for expmt. data table

    directory = directory + "kc_p_w_init_q/";
    directory = directory + sprintf("exp_%02d/", exp_num);
    file_save = directory + "kc_doublelayer_q-optim_data";

    file_save = file_save + "_" + cal_type;
    disp(directory)
    set(0,'DefaultAxesFontSize', 20); set(0,'DefaultAxesTitleFontWeight','normal');

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

    %=========================================================================
    %% main
    low_bound = 0.0;
    up_bound = 1.0;

    Tol = 1e-14;
    opts = optimset('fmincon');
    options = optimset(opts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
            'MaxFunEvals',10000,'Display','notify');

    cost_function = @(q) cost_fn(q, kc_s, kc_h);

    q0 = p;
    [q_best, min_err, exitflag, output] = fmincon(cost_function, q0, [],[],[],[],low_bound,up_bound,[],options);

    %% Generate the shapes
    pmat_total = cell(1, length(lengths));
    lengths = sort(lengths); % sort the lengths in increasing order

    % needle trajectories for all lengths
    for i = 1:length(lengths)
        Li = lengths(i);
        Ni = N_arclengths(i);

        % add the shape to a cell
        kc_s_i = kc_s * (L/Li)^p;
        kc_h_i = kc_h * (L/Li)^p;
        w_init = w_init_ref * (L/Li)^q_best;
        [~,pmati,~] = fn_intgEP_v3_2layers(w_init,kc_s_i,kc_h_i,z_crit,theta0,0,ds,Ni,B,Binv);

        % add the shape to a cell
        pmat_total{i} = pmati;

    end

    %% Calculate the deviation errors
    error_pmat_total = error_s_positions(pmat_total);
    s = 0:ds:lengths(end-1); % arclengths


    %% plot of shapes
    close all;
    figure(1);
    set(gcf,'units', 'normalized', 'position', [1/8, 1/4, 3/8, 3/8]);
    figure(2);
    set(gcf,'units', 'normalized', 'position', [1/2, 1/4, 3/8, 3/8]);

    for i = length(lengths):-1:1
        pmati = pmat_total{i};

        figure(1);
        plot3(pmati(3,:), pmati(1,:), pmati(2,:),'linewidth',2, ...
            'DisplayName', sprintf('L = %.1f mm', lengths(i))); hold on;


        figure(2);
        subplot(2,1,1);
        plot(pmati(3,:), pmati(1,:),'linewidth',2, ...
           'DisplayName', sprintf('L = %.1f mm', lengths(i))); hold on;


        subplot(2,1,2);
        plot(pmati(3,:), pmati(2,:),'linewidth',2, ...
           'DisplayName', sprintf('L = %.1f mm', lengths(i))); hold on;
        xlabel('z [mm]'); ylabel('y [mm]');

    end

    % figure titling
    f1 = figure(1);
    % patch for boundary of tissue
    [z_crit_found, idx] = min(abs(pmati(3,:)-z_crit)); % find the z_crit index
    S1 = [z_crit,        z_crit;
          pmati(1,idx) - 5, pmati(1,idx) + 5;
          pmati(1,idx) - 5 , pmati(2,idx) + 1];
    S2 = S1;
    S2(2,:) = S1(2,end:-1:1);
    S = [S1(:,1) S2(:,1) S1(:,2) S2(:,2)];

    patch(S(1,:), S(2,:), S(3,:), 'r','DisplayName', 'Tissue Boundary');

    set(gcf,'units', 'normalized', 'position', [1/10, 1/10, 6/8, 6/8]);
    hold off;
    xlabel('z [mm]'); ylabel('x [mm]'); zlabel('y [mm]')
    title(sprintf('3-D double-layer insertion: best q = %.5f', q_best))
    legend('FontSize',16);
    grid on; axis equal;

    f2 = figure(2);
    set(gcf,'units', 'normalized', 'position', [1/8, 1/8, 6/8, 6/8]);
%     sgtitle(sprintf('2-D double-layer insertion: best q = %.5f', q_best))

    subplot(2,1,1);
    ax1 = gca;
    ax1.XColor = 'k'; ax1.YColor = 'k';
    % plot the tissue boundary
    xline(z_crit, 'r--', 'Tissue Boundary','DisplayName','');
%     hold off;
    xlabel('x-z axis view, z [mm]'); ylabel('x [mm]');
%     xlabel('z [mm]'); ylabel('x [mm]');
%     title('x-z axis view');
    lgd = legend('FontSize', 16, 'location', 'nw');
    grid on; axis equal;
    
    % plot the x-axis error
    ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
    line(s, abs(error_pmat_total(1,:)), 'LineWidth', 2, 'Parent', ax2, 'Color','r');
    xlabel('s [mm]'); ylabel('Avg. Dev. [mm]');
    ylim([0 .5]);
    hold off;

    subplot(2,1,2);
    ax1 = gca;
    ax1.XColor = 'k'; ax1.YColor = 'k';
    % plot the tissue boundary
    xline(z_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
%     hold off;
    xlabel('y-z axis view, z [mm]'); ylabel('y [mm]');
%     title('y-z axis view');
    grid on; axis equal;
    
    % plot the y-axis error
    ax2 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none');
    line(s, abs(error_pmat_total(2,:)), 'LineWidth', 2, 'Parent', ax2, 'Color','r');
    xlabel('s [mm]'); ylabel('Avg. Dev. [mm]');
    ylim([0 .5]);
    hold off;
    
    return;

    %% plot of errors
    error_pmat_total_norm = vecnorm(error_pmat_total); % L2 norm values

    ferr = figure(3);
    set(gcf,'units', 'normalized', 'position', [1/10, 1/10, 6/8, 6/8]);
    subplot(2,2,1);
    plot(s, error_pmat_total_norm, 'LineWidth', 2);
    xline(z_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
    xlabel('s [mm]'); ylabel('Avg. Deviation Error [mm]');
    title("Norm of Arclength Errors");
    ylim([0 .5]);
    grid on; 

    subplot(2,2,2);
    plot(s, abs(error_pmat_total(1,:)), 'LineWidth', 2);
    xline(z_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
    ylim([0 .5]);
    title('x-axis error');
    xlabel('s [mm]');
    ylim([0 .5]);
    grid on;

    subplot(2,2,3);
    plot(s, abs(error_pmat_total(2,:)), 'LineWidth', 2);
    xline(z_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
    title('y-axis error');
    xlabel('s [mm]'); ylabel('Avg. Deviation Error [mm]');
    ylim([0 .5]);
    grid on;

    subplot(2,2,4);
    plot(s, abs(error_pmat_total(3,:)), 'LineWidth', 2);
    xline(z_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
    title('z-axis error');
    xlabel('s [mm]');
    ylim([0 .5]);
    grid on;

    sgtitle('Average Prediction Deviation Error');

    %% saving
    file_save = sprintf(file_save, p);

    msg = "";
    msg = msg + "fixed p-value and optimizing q-value\n";
    msg = msg + sprintf('kc_layer1(90 mm), %.6f\n', kc_s);
    msg = msg + sprintf('kc_layer2(90 mm), %.6f\n', kc_h);
    msg = msg + sprintf('w_init(90 mm), [ %.5f; %.5f; %.5f ]\n', w_init_ref);
    msg = msg + sprintf('p, %.5f\n', p);
    msg = msg + sprintf('best q, %.5f\n', q_best);
    msg = msg + sprintf('TipAError (mm^2), %.3f\n', min_err);

    disp(msg);

    if save_bool
        % figure
        savefig(f1, file_save + '_3D.fig');
        fprintf('Saved figure: %s\n', file_save + '_3D.fig');
        saveas(f1, file_save + '_3D.png');
        fprintf('Saved image: %s\n', file_save + '_3D.png');
        disp(' ');

        savefig(f2, file_save + '_2D.fig');
        fprintf('Saved figure: %s\n', file_save + '_2D.fig');
        saveas(f2, file_save + '_2D.png');
        fprintf('Saved image: %s\n', file_save + '_2D.png');
        disp(' ');

        savefig(ferr, file_save + '_err.fig');
        fprintf('Saved figure: %s\n', file_save + '_err.fig');
        saveas(ferr, file_save + '_err.png');
        fprintf('Saved image: %s\n', file_save + '_err.png');
        disp(' ');

        % csv file
        fid = fopen(file_save + '_err.csv', 'w');

        fprintf(fid, msg);

        fclose(fid);
        fprintf('Wrote file: %s\n', file_save + '_err.csv');

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