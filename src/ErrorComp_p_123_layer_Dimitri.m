%% ErrorComp_p_2layer_Dimitri.m

% clear all

global N_arclengths lengths B Binv ds L z1_crit z2_crit s1_crit s2_crit theta0

save_bool = true;
file_name = "Dimitri/Data/" + "kc-ideal_123layer_p-optim_error";

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

lengths = 90:15:150;
N_arclengths = floor(lengths/ds) + 1;

%% intrinsic curvature
% generalized soft, medium, hard kc values
kc_s = 0.002; kc_m = 0.003; kc_h = 0.004;

% first reference value when L = 90 mm
kc1 = kc_s;
kc2 = kc_h;
kc3 = kc1;

theta0 = 0;

%% determine s_crit (from 90 mm homogeneous single bend insertion)
% z values for the boundaries
z1_crit = 35; % in mm
z2_crit = z1_crit + 40; % in mm

% s values for the boundaries
[s1_crit, s2_crit] = determine_s_crit(kc1, kc2, z1_crit, z2_crit);

%% p values as an array
p_start = 0;
p_end = 2;
pv = p_start : 0.05 : p_end;

%=========================================================================
%% main loop
data = table('Size', [0, 7], 'VariableTypes', ["string", "double", "double", "double", "double", "double", "double"], ...
    'VariableNames', ["Layer", "p", "kc1", "kc2", "kc3", "cost", "cost_total"]);
for i = 1:length(pv)
    p = pv(i);

    [cost_total, cost_l1, cost_l2, cost_l3] = cost_fn(p, kc1, kc2, kc3);
    
    data = [data; 
            {"single", p, kc1,  -1,  -1, cost_l1, cost_total};
            {"double", p, kc1, kc2,  -1, cost_l2, cost_total};
            {"three",  p, kc1, kc2, kc3, cost_l3, cost_total}];
        
end

%% Display best p-value
[min_err, min_err_idx] = min(data.cost_total);
best_p = data.p(min_err_idx);

disp(best_p);

%% plots
f_all = figure(1);
set(f_all, 'units', 'normalized', 'position', [0.15 0.15 0.7 0.7]);
y_limit = [0, 1.1*max(data.cost)];
% single-layer
subplot(2, 3, 1);
plot(data.p(data.Layer == "single"), data.cost(data.Layer == "single"),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p'); ylabel('Error [mm^2]'); 
grid on;
title('single-layer');
ylim(y_limit);

subplot(2, 3, 2);
plot(data.p(data.Layer == "double"), data.cost(data.Layer == "double"),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p');
grid on;
title('double-layer');
ylim(y_limit);

subplot(2, 3, 3);
plot(data.p(data.Layer == "three"), data.cost(data.Layer == "three"),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p');
grid on;
title('three-layer');
ylim(y_limit);

subplot(2, 3, 5);
plot(data.p(data.Layer == "single"), data.cost_total(data.Layer == "single"),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p'); ylabel('Error [mm^2]');
grid on;
title('conjoined error');

f_cum = figure(2);
plot(data.p(data.Layer == "single"), data.cost_total(data.Layer == "single"),...
    'k.', 'markersize', 12);
xline(best_p, 'r-', 'LineWidth', 1.75);
xlabel('p'); ylabel('Error [mm^2]');
grid on;
title('conjoined error');

%% Save the Results
if save_bool
    writetable(data, file_name + '.csv')
    fprintf('Wrote csv file: %s\n', file_name + '.csv');

    % save just the conjoined error
    savefig(f_cum, file_name + '_cum.fig');
    fprintf('Wrote Figure: %s\n', file_name + '_cum.png');

    saveas(f_cum, file_name + '_cum.png');
	fprintf('Wrote Image: %s\n', file_name + '_cum.png');

    % save all plots error plot
    savefig(f_all, file_name + '_all.fig');
    fprintf('Wrote Figure: %s\n', file_name +  '_all.png');

    saveas(f_all, file_name + '_all.png');
	fprintf('Wrote Image: %s\n', file_name +  '_all.png');
    
end

%% Functions
% cost function for the 1 & 2-layer cases
function [cost, cost_1layer, cost_2layer, cost_3layer] = cost_fn(p, kc1, kc2, kc3)

    cost_1layer = cost_fn_1layer(p, kc1);
    cost_2layer = cost_fn_2layer(p, kc1, kc2);
    cost_3layer = cost_fn_3layer(p, kc1, kc2, kc3);
    cost = cost_1layer + cost_2layer + cost_3layer;

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
    global L lengths ds theta0 N_arclengths B Binv z1_crit s1_crit
    
    pmat_total = cell(1,length(lengths));
    
    % needle trajectories for all lengths
    for i = 1:length(lengths)
        Li = lengths(i);
        Ni = N_arclengths(i);
        
        % generate the shape
        kc1_i = kc1 * (L/Li)^p;
        kc2_i = kc2 * (L/Li)^p;
        k0 = kc1_i*(s1_crit/Li)^2 + kc2_i*(1 - s1_crit^2/Li^2);
        w_init = [k0;0;0]; 
        [~,pmati,~] = fn_intgEP_v3_2layers(w_init,kc1_i,kc2_i,z1_crit,theta0,0,ds,Ni,B,Binv);
        
        % add the shape to a cell
        pmat_total{i} = pmati;
        
    end
    
    cost = TipAerror_Dimitri(pmat_total, lengths);
    
end

% cost function for optimization
function cost = cost_fn_3layer(p, kc1, kc2, kc3)    
    global L lengths ds theta0 N_arclengths B Binv s1_crit s2_crit
    
    pmat_total = cell(1,length(lengths));
    
    % needle trajectories for all lengths
    for i = 1:length(lengths)
        Li = lengths(i);
        Ni = N_arclengths(i);
        
        % generate the shape
        kc1_i = kc1 * (L/Li)^p;
        kc2_i = kc2 * (L/Li)^p;
        kc3_i = kc3 * (L/Li)^p;
%         w_init = w_init_ref * (L/Li)^p;
        w_init = (kc1_i*(s1_crit).^2/L^2 + kc2_i*(s2_crit - s1_crit)/L*(s2_crit + s1_crit)/L ...
                + kc3_i*(1 - s2_crit/L)*(1 + s2_crit/L)) * [1; 0; 0]; % ideal
        [~,pmati,~] = fn_intgEP_3layers_s_crit_Dimitri(w_init, kc1_i, kc2_i, kc3_i, ...
                     theta0, 0, ds, s1_crit, s2_crit, Li, B, Binv);
        
        % add the shape to a cell
        pmat_total{i} = pmati;
        
    end
    
    cost = TipAerror_Dimitri(pmat_total, lengths);
    
end

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