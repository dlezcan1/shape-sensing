%% ErrorComp_p_2layer_Dimitri.m

% clear all

global N_arclengths lengths B Binv ds L z_crit s_crit theta0

save_bool = true;

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

z_crit = 35;

lengths = [60, 90, 120, 150, 180];
N_arclengths = floor(lengths/ds) + 1;
theta0 = 0;

%% intrinsic curvature
% use this for two-layer transition
kc1 = 0.002;
kc2 = 0.004;

%% determine s_crit (from 90 mm homogeneous single bend insertion)
[~,pos1] = fn_intgEP_v1_1layer([kc1;0;0],kc1,0,0,ds,N,B,Binv);
ix_crit_v = find(abs(pos1(3,:) - z_crit) <= ds/2); 
ix_crit = ix_crit_v(1);
s_crit = s(ix_crit);

%% p values as an array
p_start = 0.4;
p_end = 0.7;
Np = 51;
dp = (p_end-p_start)/(Np-1);
pv = p_start:dp:p_end;

%=========================================================================
%% main loop
data = table('Size', [0, 6], 'VariableTypes', ["string", "double", "double", "double", "double", "double"], ...
    'VariableNames', ["Layer", "p", "kc1", "kc2", "cost", "cost_total"]);
for i = 1:Np
    p = pv(i);

    [cost_total, cost_l1_soft, cost_l1_hard, cost_l2, cost_l2_rev] = cost_fn(p, kc1, kc2);
    
    data = [data; 
            {"single", p, kc1, -1, cost_l1_soft, cost_total};
            {"single", p, kc2, -1, cost_l1_hard, cost_total};
            {"double", p, kc1, kc2, cost_l2, cost_total};
            {"double", p, kc2, kc1, cost_l2_rev, cost_total}];
        
end

%% Display best p-value
[min_err, min_err_idx] = min(data.cost_total);
best_p = data.p(min_err_idx);

disp(best_p);

%% plots
f_all = figure(1);

% single-layer
% kc1
subplot(2, 2, 1);
plot(data.p(data.Layer == "single" & data.kc1 == kc1), data.cost(data.Layer == "single" & data.kc1 == kc1),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p'); ylabel('Error [mm^2]'); 
grid on;
title('single-layer: soft');

% kc2
subplot(2, 2, 2);
plot(data.p(data.Layer == "single" & data.kc1 == kc2), data.cost(data.Layer == "single" & data.kc1 == kc2),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p'); ylabel('Error [mm^2]'); 
grid on;
title('single-layer: hard');

% double-layer
% kc1 - kc2
subplot(2, 2, 3);
plot(data.p(data.Layer == "double" & data.kc1 == kc1), data.cost(data.Layer == "double" & data.kc1 == kc1),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p');
grid on;
title('double-layer: soft-to-hard');

% kc2 - kc1
subplot(2, 2, 4);
plot(data.p(data.Layer == "double" & data.kc1 == kc2), data.cost(data.Layer == "double" & data.kc1 == kc2),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p');
grid on;
title('double-layer: hard-to-soft');

% conjoined error
% subplot(2, 2, 5);
% plot(data.p(data.Layer == "single" & data.kc1 == kc1), data.cost_total(data.Layer == "single" & data.kc1 == kc1),...
%     'k.', 'markersize', 12);
% xline(best_p, 'r-', 'LineWidth', 1.75);
% xlabel('p');
% grid on;
% title('conjoined error');

f_sum = figure(2);
plot(data.p(data.Layer == "single"), data.cost_total(data.Layer == "single"),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p'); ylabel('Error [mm^2]');
grid on;
title('conjoined error');

%% Save the Results
file_name = "Dimitri/Data/" + "kc-ideal_p-optim_error";

if save_bool
    writetable(data, file_name + '.csv')
    fprintf('Wrote csv file: %s\n', file_name + '.csv');

    savefig(f_sum, file_name + '.fig');
    fprintf('Wrote Figure: %s\n', file_name + '.fig');

    saveas(f_sum, file_name + '.png');
    fprintf('Wrote Image: %s\n', file_name + '.png');
    
    savefig(f_all, file_name + '_all.fig');
    fprintf('Wrote Figure: %s\n', file_name + '_all.fig');

    saveas(f_all, file_name + '_all.png');
    fprintf('Wrote Image: %s\n', file_name + '_all.png');
    
    
end

%% Functions
% cost function for the 1 & 2-layer cases
function [cost, cost_1layer_soft, cost_1layer_hard, cost_2layer, cost_2layer_rev] = cost_fn(p, kc1, kc2)

    cost_1layer_soft = cost_fn_1layer(p, kc1);
    cost_1layer_hard = cost_fn_1layer(p, kc2);
    cost_2layer = cost_fn_2layer(p, kc1, kc2);
    cost_2layer_rev = cost_fn_2layer(p, kc2, kc1); % switch the layers
    cost = cost_1layer_soft + cost_1layer_hard + cost_2layer + cost_2layer_rev;

end

function cost = cost_fn_1layer(p, kc)
    global lengths L ds B Binv N_arclengths
    
    pmat_total = cell(1,length(lengths));
    
    % generate the positions for ideal insertion
    for i = 1:length(lengths)
       Li = lengths(i);
       Ni = N_arclengths(i);
       
       kc_i = kc * (L/Li)^p;
       w_init = [kc_i;0;0];
       [~, pmat_i, ~] = fn_intgEP_v1_1layer(w_init,kc_i,0,0,ds,Ni,B,Binv);
       
       pmat_total{i} = pmat_i;
       
    end
    
    cost = TipAerror_Dimitri(pmat_total, lengths);
    
end


function cost = cost_fn_2layer(p, kc1, kc2)
    global L lengths ds theta0 N_arclengths B Binv z_crit s_crit
    
    pmat_total = cell(1,length(lengths));
    
    % needle trajectories for all lengths
    for i = 1:length(lengths)
        Li = lengths(i);
        Ni = N_arclengths(i);
        
        % generate the shape
        kc1_i = kc1 * (L/Li)^p;
        kc2_i = kc2 * (L/Li)^p;
        k0 = kc1_i*(s_crit/Li)^2 + kc2_i*(1 - s_crit^2/Li^2);
        w_init = [k0;0;0]; 
        [~,pmati,~] = fn_intgEP_v3_2layers(w_init,kc1_i,kc2_i,z_crit,theta0,0,ds,Ni,B,Binv);
        
        % add the shape to a cell
        pmat_total{i} = pmati;
        
    end
    
    cost = TipAerror_Dimitri(pmat_total, lengths);
    
end