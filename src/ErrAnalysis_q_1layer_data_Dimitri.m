%% ErrAnalysis_q_1layer_data_Dimitri.m
%
% from ideal_singlebendinsertion_2020_0413.m
%
% used for the simulation of data-based single bend needle insertion.
%
% - written by Jin Seob Kim
% - edited  by Dimitri Lezcano

global p L L1 L2 L3 L4 L5 ds theta0 N1 N2 N3 N4 N5 B Binv

%% saving set-up
save_bool = true;
directory = "Dimitri/Data/";
directory = directory + "Archive/1-Layer/TipAError_Cost_Ideal/";
file_save = directory + "kc-q_singlelayer_p-%.3f_ideal";

disp(directory)

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
theta0 = 0.0;
kc = 0.003;

%% generate arclength array
ds = 0.5; % in mm
L1 = 90; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 105; N2 = L2/ds+1; s2 = linspace(0,L2,N2);]
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

lengths = 90:15:150;
N_arclengths = floor(lengths/ds) + 1;

%% dL and n values as an array
p = 0.6;
dL_v = .5:.5:15; N_dL = length(dL_v);
n_v = 2:10;  N_n = length(n_v);

%=========================================================================
%% main loop
data_tbl = array2table(zeros(length(dL_v)*length(n_v), 5), ...
    'VariableNames', {'dL', 'n', 'q', 'TipAerror_mm2', 'Time_s'});

h = waitbar(0, 'Please wait...', 'Name', 'dL, n | Progress');
for i = 1:N_dL
    for j = 1:N_n
        % parameters
        dL = dL_v(i);
        n = n_v(j);
        
        % begin optimization
        tic; 
        low_bound = 0.0;
        up_bound = 1.0;

        Tol = 1e-14;
        opts = optimset('fmincon');
        options = optimset(opts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
                'MaxFunEvals',10000,'Display','notify');

        cost_function = @(q) cost_fn(q, kc, dL, n);

        q0 = p;
        [q_best, fval, exitflag, output] = fmincon(cost_function, q0, [],[],[],[],low_bound,up_bound,[],options);
        t = toc;
        
        % add the data to a table
        idx = sub2ind([ N_n, N_dL ], j, i);
        data_tbl(idx,:) = {dL, n, q_best, fval, t};


        waitbar(idx/(N_dL * N_n), h, sprintf("%.2f%% complete", 100*idx/(N_dL * N_n)));
    
    end
    
end
close(h);

%% plots
% grid the data
grid_dL = reshape(data_tbl.dL, N_dL, []);
grid_n = reshape(data_tbl.n, N_dL, []);
grid_q = reshape(data_tbl.q, N_dL, []);
grid_err = reshape(data_tbl.TipAerror_mm2, N_dL, []);
grid_time = reshape(data_tbl.Time_s, N_dL, []);

[grid_int_dL, grid_int_n] = meshgrid(data_tbl.dL(1):.1:data_tbl.dL(end), data_tbl.n(1):.1:data_tbl.n(end));
grid_int_q = griddata(grid_dL, grid_n, grid_q, grid_int_dL, grid_int_n);
grid_int_err = griddata(grid_dL, grid_n, grid_err, grid_int_dL, grid_int_n);
grid_int_time = griddata(grid_dL, grid_n, grid_time, grid_int_dL, grid_int_n);

grid_int_err_per_length = grid_int_err ./ (grid_int_dL .* grid_int_n + L);
grid_int_err_per_length_per_n = grid_int_err_per_length./ grid_int_n;

f1 = figure(1);
set(gcf,'units', 'normalized', 'position', [0, 0, 1, 1]);
subplot(2,2,1);
surf(grid_int_dL, grid_int_n, grid_int_err, 'edgecolor', 'none');
xlabel('dL [mm]'); ylabel('# of insertions'); zlabel('Area Deviation Error [mm^2]')
title('Area Deviation Error');
grid on

f1 = figure(1);
subplot(2,2,2);
surf(grid_int_dL, grid_int_n, grid_int_err_per_length, 'edgecolor', 'none');
xlabel('dL [mm]'); ylabel('# of insertions'); zlabel('Avg. Deviation [mm]')
title('Average Deviation');
grid on

subplot(2,2,3);
surf(grid_int_dL, grid_int_n, grid_int_time, 'edgecolor', 'none');
xlabel('dL [mm]'); ylabel('# of insertions'); zlabel('Run time [s]')
title('Time for Optimization');
grid on

subplot(2,2,4);
surf(grid_int_dL, grid_int_n, grid_int_err_per_length_per_n, 'edgecolor', 'none');
xlabel('dL [mm]'); ylabel('# of insertions'); zlabel('Avg. Deviation/Insertion [mm]')
title('Avg. Deviation per Simulated Insertion [mm]');
grid on

sgtitle(sprintf("Ideal Single-Layer Insertion: p = %.4f", p));

%% saving
file_save = sprintf(file_save, p);

if save_bool
    % figure
    savefig(f1, file_save + '_dL-n_optim.fig');
    fprintf('Saved figure: %s\n', file_save + '_dL-n_optim.fig');
    saveas(f1, file_save + '_dL-n_optim.png');
    fprintf('Saved image: %s\n', file_save + '_dL-n_optim.png');
    
    % data table
    writetable(data_tbl, file_save + '_dL-n_optim.csv')
    fprintf('Saved data table: %s\n', file_save + '_dL-n_optim.csv');
    disp(' ');
    
end   

%% Functions
% cost function for optimization
function cost = cost_fn(q, kc, dL, n)    
    global p L ds theta0 B Binv
    pmat_total = cell(1, n);
    lengths = zeros(1, n);
    for i = 1:n
        Li = L + dL*(i-1);
        lengths(i) = Li;
        
        % needle trajectories
        kci = kc*(L/Li)^p;
        w_init = [kc;0;0]*(L/Li)^q;
        Ni = floor(Li/ds);
        [~,pmati,~] = fn_intgEP_v1_1layer(w_init,kci,theta0,0,ds,Ni,B,Binv);
        pmat_total{i} = pmati;
    
    end

    cost = TipAerror_Dimitri(pmat_total, lengths);
    
end


