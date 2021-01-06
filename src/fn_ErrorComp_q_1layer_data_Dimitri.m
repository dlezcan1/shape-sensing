%% fn_ErrorComp_q_1layer_data_Dimitri.m
%
% from ideal_singlebendinsertion_2020_0413.m
%
% used for the simulation of data-based single bend needle insertion.
%
% - written by Jin Seob Kim
% - edited  by Dimitri Lezcano

function [pmat_total, q_best, min_err] = fn_ErrorComp_q_1layer_data_Dimitri(exp_num, cal_type, save_bool, exp_type)
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
    file_save = directory + "kc_singlelayer_q-optim_data";
    
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
    p = 0.5920;

    % data
    w_init_ref = reshape(data_row.w_init_ref,3,1);
    theta0 = 1.0*pi/180;

    %% generate arclength array
    ds = 0.5; % in mm

    lengths = 90:15:150;
    N_arclengths = lengths/ds + 1;
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
    [q_best, min_err, exitflag, output] = fmincon(cost_function, q0, [],[],[],[],low_bound,up_bound,[],options);

    %% Generate q_best insertions
    % needle trajectories
    pmat_total = cell(1, length(lengths));
    lengths = sort(lengths); % sort the lengths in increasing order
    
    % needle trajectories for all lengths
    for i = 1:length(lengths)
        Li = lengths(i);
        Ni = N_arclengths(i);
        
        % add the shape to a cell
        kci = kc_s*(L/Li)^p;
        w_init = w_init_ref * (L/Li)^q_best;
        [~,pmati,~] = fn_intgEP_v1_1layer(w_init,kci,theta0,0,ds,Ni,B,Binv);
        
        % add the shape to a cell
        pmat_total{i} = pmati;
        
    end

    %% Calculate the deviation errors
    error_pmat_total = error_s_positions(pmat_total);
    s = 0:ds:lengths(end-1); % arclengths
    
    %% plots
    close all;
    f1 = figure(1);
    set(gcf,'units', 'normalized', 'position', [1/8, 1/4, 3/8, 3/8]);
    f2 = figure(2);
    set(gcf,'units', 'normalized', 'position', [1/2, 1/4, 3/8, 3/8]);
    for i = length(pmat_total):-1:1
        figure(1);
        pmati = pmat_total{i};
        plot3(pmati(3,:), pmati(1,:), pmati(2,:), 'linewidth', 2, 'displayname', sprintf('L = %.1f', lengths(i)));
        hold on;

        figure(2);
        subplot(2,1,1);
        plot(pmati(3,:), pmati(1,:), 'linewidth', 2, 'displayname', sprintf('L = %.1f mm', lengths(i))); 
        hold on;

        subplot(2,1,2);
        plot(pmati(3,:), pmati(2,:), 'linewidth', 2, 'displayname', sprintf('L = %.1f mm', lengths(i))); 
        hold on;

    end
    figure(1);
    set(gcf,'units', 'normalized', 'position', [1/10, 1/10, 6/8, 6/8]);
    hold off;
    xlabel('z [mm]'); ylabel('x [mm]'); zlabel('y [mm]')
    title(sprintf("3-D single-layer insertion: best p' = %.5f", q_best))
    legend();
    grid on; axis equal;

    figure(2);
    set(gcf,'units', 'normalized', 'position', [1/8, 1/8, 6/8, 6/8]);
%     sgtitle(sprintf("2-D single-layer insertion: best p' = %.5f", q_best))
    subplot(2,1,1);
    ax1 = gca;
    hold off;
%     xlabel('z [mm]'); ylabel('x [mm]');
%     title('x-z axis view');
    xlabel('x-z axis view, z [mm]'); ylabel('x [mm]');
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
    
    %% plot of errors
    error_pmat_total_norm = vecnorm(error_pmat_total); % L2 norm values

    ferr = figure(3);
    set(gcf,'units', 'normalized', 'position', [1/7, 1/7, 6/8, 6/8]);
    subplot(2,2,1);
    plot(s, error_pmat_total_norm, 'LineWidth', 2);
    xlabel('s [mm]'); ylabel('Avg. Deviation Error [mm]');
    title("Norm of Arclength Errors");
    ylim([0 .5]);
    grid on; 

    subplot(2,2,2);
    plot(s, abs(error_pmat_total(1,:)), 'LineWidth', 2);
    ylim([0 .5]);
    title('x-axis error');
    xlabel('s [mm]');
    ylim([0 .5]);
    grid on;

    subplot(2,2,3);
    plot(s, abs(error_pmat_total(2,:)), 'LineWidth', 2);
    title('y-axis error');
    xlabel('s [mm]'); ylabel('Avg. Deviation Error [mm]');
    ylim([0 .5]);
    grid on;

    subplot(2,2,4);
    plot(s, abs(error_pmat_total(3,:)), 'LineWidth', 2);
    title('z-axis error');
    xlabel('s [mm]');
    ylim([0 .5]);
    grid on;

    sgtitle('Average Prediction Deviation Error');

    %% saving
    file_save = sprintf(file_save, p);

    msg = "";
    msg = msg + sprintf("fixed p-value and optimizing q-value\n");
    msg = msg + sprintf('kc(90 mm), %.6f\n', kc_s);
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
        disp(" ");

    end   
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