%% kc_q_1layer_data_Dimitri.m
%
% from ideal_singlebendinsertion_2020_0413.m
%
% used for the simulation of data-based single bend needle insertion.
%
% - written by Jin Seob Kim
% - edited  by Dimitri Lezcano

clear all;

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

%% generate arclength array
ds = 0.5; % in mm
L1 = 90; N1 = L1/ds+1; s1 = linspace(0,L1,N1);
L2 = 105; N2 = L2/ds+1; s2 = linspace(0,L2,N2);
L3 = 120; N3 = L3/ds+1; s3 = linspace(0,L3,N3);
L4 = 150; N4 = L4/ds+1; s4 = linspace(0,L4,N4);
L5 = 180; N5 = L5/ds+1; s5 = linspace(0,L5,N5);

lengths = [L1, L2, L3, L4, L5];

%% q values as an array
p = 0.6;
kc_start = 0.0015;
kc_end = 0.006;
kc_v = linspace(kc_start, kc_end, 101);

%=========================================================================
%% main loop
data_tbl = array2table(zeros(length(kc_v), 3),...
    'VariableNames', {'kc', 'q', 'TipAerror'});
h = waitbar(0, 'Please wait...');
for i = 1:length(kc_v)
    kc = kc_v(i);
    
    low_bound = 0.0;
    up_bound = 1.0;

    Tol = 1e-14;
    opts = optimset('fmincon');
    options = optimset(opts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
            'MaxFunEvals',10000,'Display','notify');
        
    cost_function = @(q) cost_fn(q, kc);
        
    q0 = p;
    [q_best, fval, exitflag, output] = fmincon(cost_function, q0, [],[],[],[],low_bound,up_bound,[],options);
    data_tbl(i,:) = {kc, q_best, fval};
    
    waitbar(i/length(kc_v), h, sprintf("%.2f%% complete", i/length(kc_v)*100));
end
close(h);


%% plots
f1 = figure(1);
set(gcf,'units', 'normalized', 'position', [1/4, 1/4, 1/2, 1/2]);
subplot(1,2,1);
plot(data_tbl.kc,data_tbl.q,'ko-','linewidth',2)
xlabel('\kappa_c')
ylabel('optimized q-value')
grid on

subplot(1,2,2);
plot(data_tbl.kc,data_tbl.TipAerror,'ko-','linewidth',2)
xlabel('\kappa_c')
ylabel('error [mm^2]')
grid on

sgtitle("ideal single-layer insertion: \kappa_c vs. q")

%% saving
file_save = sprintf(file_save, p);

if save_bool
    % figure
    savefig(f1, file_save + '.fig');
    fprintf('Saved figure: %s\n', file_save + '.fig');
    saveas(f1, file_save + '.png');
    fprintf('Saved image: %s\n', file_save + '.png');
    
    % data table
    writetable(data_tbl, file_save + '.csv')
    fprintf('Saved data table: %s\n', file_save + '.csv');
    disp(' ');
    
end   


%% Functions

% cost function for optimization
function cost = cost_fn(q, kc)    
    global p L L1 L2 L3 L4 L5 ds theta0 N1 N2 N3 N4 N5 B Binv
    
    % needle trajectories
    kc1 = kc*(L/L1)^p;
    w_init = [kc;0;0]*(L/L1)^q;
    [wv1,pmat1,Rmat1] = fn_intgEP_v1_1layer(w_init,kc1,theta0,0,ds,N1,B,Binv);

    kc2 = kc*(L/L2)^p;
    w_init = [kc;0;0]*(L/L2)^q;
    [wv2,pmat2,Rmat2] = fn_intgEP_v1_1layer(w_init,kc2,theta0,0,ds,N2,B,Binv);

    kc3 = kc*(L/L3)^p;
    w_init = [kc;0;0]*(L/L3)^q;
    [wv3,pmat3,Rmat3] = fn_intgEP_v1_1layer(w_init,kc3,theta0,0,ds,N3,B,Binv);

    kc4 = kc*(L/L4)^p;
    w_init = [kc;0;0]*(L/L4)^q;
    [wv4,pmat4,Rmat4] = fn_intgEP_v1_1layer(w_init,kc4,theta0,0,ds,N4,B,Binv);

    kc5 = kc*(L/L5)^p;
    w_init = [kc;0;0]*(L/L5)^q;
    [wv5,pmat5,Rmat5] = fn_intgEP_v1_1layer(w_init,kc5,theta0,0,ds,N5,B,Binv);
    
    % error
    pmat_total = cell(1,5);
    pmat_total{1} = pmat1;
    pmat_total{2} = pmat2;
    pmat_total{3} = pmat3;
    pmat_total{4} = pmat4;
    pmat_total{5} = pmat5;
    
    cost = TipAerror(pmat_total);
    
end


