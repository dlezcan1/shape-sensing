%% ErrorComp_p_3layer_Dimitri.m
%
% used for the simulation of data-based single bend needle insertion.
%
% - written by Dimitri Lezcano

global L lengths N ds theta0 B Binv w_init_ref s1_crit s2_crit N_arclengths z1_crit z2_crit

% saving set-up
save_bool = false;
directory = "Dimitri/Data/";
directory = directory + "Archive/3-Layer/TipAError_Cost_Ideal/";
file_save = directory + "kc_3layer_p-err_ideal";

if save_bool
    fprintf("Save Directory: %s\n", directory)
    disp(" ");
end

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
% generalized soft, medium, hard kc values
kc_s = 0.002; kc_m = 0.003; kc_h = 0.004;

% first reference value when L = 90 mm
kc1 = kc_s;
kc2 = kc_h;
kc3 = kc1;

theta0 = 0;

%% generate arclength array
lengths = 90:15:150;
N_arclengths = floor(lengths/ds) + 1;

%% determine s_crit (from 90 mm homogeneous single bend insertion)
% z values for the boundaries
z1_crit = 35; % in mm
z2_crit = z1_crit + 40; % in mm

% s values for the boundaries
[s1_crit, s2_crit] = determine_s_crit(kc1, kc2, z1_crit, z2_crit);

%% Shape-sensing initial conditions
% ideal
w_init_ref = (kc1*(s1_crit).^2/L^2 + kc2*(s2_crit - s1_crit)/L*(s2_crit + s1_crit)/L ...
        + kc3*(1 - s2_crit/L).*(1 + s2_crit/L)) * [1; 0; 0]; 

%% Values to run over optimization (p)
p_start = 0.4;
p_end = 0.7;
Np = 51;
dp = (p_end-p_start)/(Np-1);
pv = 0:0.05:2;

%=========================================================================
%% Main error
cost_func = @(p1, p2) cost_fn(p1, p2, kc1, kc2, kc3); % cost function handle
data_tbl = table('size', [0, 7], 'variabletypes', ["string", "double", "double", "double", "double", "double", "double"], ...
    'variablenames', ["Layers", "p12", "p3", "kc1", "kc2", "kc3", "cost"]);

for i = 1:length(pv)
    for j = 1:length(pv)
    
    % get the p-value
    p_i = pv(i);
    p_j = pv(j);
    
    % compute the cost
    cost = cost_func(p_i, p_j);
    
    data_tbl = [data_tbl;
                {"three", p_i, p_j, kc1, kc2, kc3, cost}];
    
    end
end

%% determine the best values
[~, min_idx] = min(data_tbl.cost);
p_best = data_tbl.p(min_idx);
fprintf('p_bstg = %f\n', p_best);
% p1_best = data_tbl.p1(min_idx);
% p23_best = data_tbl.p23(min_idx);
% fprintf("p1_best: %f\np23_best: %f\n", p1_best, p23_best);

%% Generate the shapes
pmat_total = cell(1, length(lengths));
lengths = sort(lengths); % sort the lengths in increasing order

% needle trajectories for all lengths
for i = 1:length(lengths)
    Li = lengths(i);
    Ni = N_arclengths(i);

    % generate the shape
    kc1_i = kc1 * (L/Li)^p1_best;
    kc2_i = kc2 * (L/Li)^p1_best;
    kc3_i = kc3 * (L/Li)^p23_best;
%     w_init = w_init_ref * (L/Li)^p_best;
    w_init = kc1_i*(s1_crit).^2/L^2 + kc2_i*(s2_crit - s1_crit)/L*(s2_crit + s1_crit)/L ...
        + kc3_i*(1 - s2_crit/L)*(1 + s2_crit/L) * [1;0;0];
    [~,pmati,~] = fn_intgEP_3layers_s_crit_Dimitri(w_init, kc1_i, kc2_i, kc3_i, ...
                     theta0, 0, ds, s1_crit, s2_crit, Li, B, Binv);

    % add the shape to a cell
    pmat_total{i} = pmati;

end

%% Calculate the errors
% calculate the deviation error
error_pmat_total = error_s_positions(pmat_total);

% calculate the TipAerror
tipAerr = TipAerror_Dimitri(pmat_total, lengths);
    
%% plot of shapes
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
    title('y-z axis view');
    
end

% figure titling
f1 = figure(1);
% patches for tissue boudaries
[~, idx1] = min(abs(pmati(3,:)-z1_crit)); % find the z_crit index
S1 = [z1_crit,        z1_crit;
      pmati(1,idx1) - 5, pmati(1,idx1) + 5;
      pmati(1,idx1) - 5 , pmati(2,idx1) + 1];
S2 = S1;
S2(2,:) = S1(2,end:-1:1);
S_l1 = [S1(:,1) S2(:,1) S1(:,2) S2(:,2)];

[~, idx2] = min(abs(pmati(3,:)-z2_crit)); % find the z_crit index
S1 = [z2_crit,        z2_crit;
      pmati(1,idx2) - 5, pmati(1,idx2) + 5;
      pmati(1,idx2) - 5 , pmati(2,idx2) + 1];
S2 = S1;
S2(2,:) = S1(2,end:-1:1);
S_l2 = [S1(:,1) S2(:,1) S1(:,2) S2(:,2)];

% plot the tissue boudary patches
patch(S_l1(1,:), S_l1(2,:), S_l1(3,:), 'r','DisplayName', 'Tissue Boundary');
patch(S_l2(1,:), S_l2(2,:), S_l2(3,:), 'r','DisplayName', 'Tissue Boundary');

set(gcf,'units', 'normalized', 'position', [1/10, 1/10, 6/8, 6/8]);
hold off;
xlabel('z [mm]'); ylabel('x [mm]'); zlabel('y [mm]')
title(sprintf('3-D triple-layer insertion: p12 = %.5f, p3 = %.5f', p1_best, p23_best))
legend('FontSize',16);
grid on; axis equal;

f2 = figure(2);
set(gcf,'units', 'normalized', 'position', [1/8, 1/8, 6/8, 6/8]);
sgtitle(sprintf('2-D triple-layer insertion: p12 = %.5f, p3 = %.5f', p1_best, p23_best))

subplot(2,1,1);
% plot the tissue boundary
xline(z1_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
xline(z2_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
hold off;
xlabel('z [mm]'); ylabel('x [mm]');
title('x-z axis view');
legend('FontSize', 16,'location','northeastoutside');
grid on; axis equal;

subplot(2,1,2);
% plot the tissue boundary
xline(z1_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
xline(z2_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
hold off;
xlabel('z [mm]'); ylabel('y [mm]');
title('y-z axis view');
grid on; axis equal;

%% plot of errors
error_pmat_total_norm = vecnorm(error_pmat_total); % L2 norm values
s = 0:ds:lengths(end-1); % arclengths

ferr = figure(3);
sgtitle('Average Prediction Deviation Error');
set(gcf,'units', 'normalized', 'position', [1/10, 1/10, 6/8, 6/8]);

err_lim = ceil(10*max([0.5, error_pmat_total_norm]))/10;

% norm of errors
subplot(2,2,1);
plot(s, error_pmat_total_norm, 'LineWidth', 2);

xline(z1_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
xline(z2_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');

xlabel('s [mm]'); ylabel('Avg. Deviation Error [mm]');
title("Norm of Arclength Errors");
ylim([0 err_lim]);
grid on; 

% x-axis errors
subplot(2,2,2);
plot(s, abs(error_pmat_total(1,:)), 'LineWidth', 2);

xline(z1_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
xline(z2_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');

ylim([0 .5]);
title('x-axis error');
xlabel('s [mm]');
ylim([0 err_lim]);
grid on;

% y-axis errors
subplot(2,2,3);
plot(s, abs(error_pmat_total(2,:)), 'LineWidth', 2);

xline(z1_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
xline(z2_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');

title('y-axis error');
xlabel('s [mm]'); ylabel('Avg. Deviation Error [mm]');
ylim([0 err_lim]);
grid on;

% z-axis errors
subplot(2,2,4);
plot(s, abs(error_pmat_total(3,:)), 'LineWidth', 2);

xline(z1_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');
xline(z2_crit, 'r--', 'Tissue Boundary','DisplayName','Tissue Boundary');

title('z-axis error');
xlabel('s [mm]');
ylim([0 err_lim]);
grid on;

%% plot the p vs. error values
f_perr = figure(4);

% [p1_grid, p23_grid] = meshgrid(unique(data_tbl.p1), unique(data_tbl.p23));
% cost_grid = reshape(data_tbl.cost, size(p1_grid, 1), size(p23_grid, 2))';
% 
% surf(p1_grid, p23_grid, cost_grid);
% xlabel('p1'); ylabel('p23'); zlabel('Error [mm^2]');
% view([0, 90])


plot(data_tbl.p, data_tbl.cost, 'k.', 'markersize', 12);
xline(p_best, 'r-', 'LineWidth', 1.75);
xlabel('p'); ylabel('Error [mm^2]');

grid on;
title('three-layer');


%% saving
file_save = sprintf(file_save, p);

msg = "";
msg = msg + sprintf("Ideal Insertion\n");
msg = msg + sprintf('kc_layer1(90 mm), %.6f\n', kc1);
msg = msg + sprintf('kc_layer2(90 mm), %.6f\n', kc2);
msg = msg + sprintf('kc_layer3(90 mm), %.6f\n', kc3);
% msg = msg + sprintf('w_init(90 mm), [ %.5f; %.5f; %.5f ]\n', w_init_ref);
msg = msg + sprintf('p_best, %.5f\n', p1_best);
% msg = msg + sprintf('best q, %.5f\n', q_best);
msg = msg + sprintf('TipAError (mm^2), %.3f\n', tipAerr);
msg = msg + sprintf('Avg. Deviation (mm), %.5f\n', tipAerr/lengths(end-1));

fprintf("Save filename base: " + file_save + "\n" + msg);
disp(" ");

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
    
    % write table
    writetable(data_tbl, file_save + '_p-err.csv');
    fprintf('Wrote table to file: %s\n', file_save + '_p-err.csv');

end

%% Functions
% function to generate the s_crit values
function [s1_crit, s2_crit] = determine_s_crit(kc1, kc2, z1_crit, z2_crit)
    global theta0 ds B Binv N L
    
    s = 0:ds:L;
    % calculate ideal w_init 
    % determine 
    w_init_1ideal = kc1 * [1; 0; 0];
    [~,p1] = fn_intgEP_v1_1layer(w_init_1ideal,kc1,theta0,0,ds,N,B,Binv);
    ix_crit_v = find(abs(p1(3,:) - z1_crit) <= ds/2); 
    ix_crit = ix_crit_v(1);
    s_crit = s(ix_crit);

    w_init_2ideal = (kc1*(s_crit).^2/L^2 + kc2*(s_crit/L)*(1 + s_crit/L)).*[1; 0; 0];    
    
    % determine s1_crit and s2_crit (90 mm)
    [~, p] = fn_intgEP_v1_2layers(w_init_2ideal, kc1, kc2, z1_crit, theta0, 0, ds, N, B, Binv);

    ix_crit_v1 = find(abs(p(3,:) - z1_crit) <= ds/2); 
    ix_crit1 = ix_crit_v1(1);
    s1_crit = s(ix_crit1); % return value

    ix_crit_v2 = find(abs(p(3,:) - z2_crit) <= ds/2); 
    ix_crit2 = ix_crit_v2(1);
    s2_crit = s(ix_crit2); % return value

end

% cost function for optimization
function cost = cost_fn(p12, p3, kc1, kc2, kc3)    
    global L lengths ds theta0 N_arclengths B Binv s1_crit s2_crit
    
    pmat_total = cell(1,length(lengths));
    
    % needle trajectories for all lengths
    for i = 1:length(lengths)
        Li = lengths(i);
        Ni = N_arclengths(i);
        
        % generate the shape
        kc1_i = kc1 * (L/Li)^p12;
        kc2_i = kc2 * (L/Li)^p12;
        kc3_i = kc3 * (L/Li)^p3;
%         w_init = w_init_ref * (L/Li)^p;
        w_init = (kc1_i*(s1_crit).^2/L^2 + kc2_i*(s2_crit - s1_crit)/L*(s2_crit + s1_crit)/L ...
                + kc3_i*(1 - s2_crit/L)*(1 + s2_crit/L)) * [1; 0; 0]; % ideal
        [~,pmati,~] = fn_intgEP_3layers_s_crit_Dimitri(w_init, kc1_i, kc2_i, kc3_i, ...
                     theta0, 0, ds, s1_crit, s2_crit, Li, B, Binv);
        
        % add the shape to a cell
        pmat_total{i} = pmati;
        
    end
    
    cost = TipAerror_ref_Dimitri(pmat_total, lengths, L);
    
end