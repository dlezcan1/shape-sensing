%% Divide_Conquer_2to1_layers.m
%
% This script is to test a method used to minimize the error between kappa_0
%  curves from 2 to 1-layers
%
% - written by: Dimitri Lezcano


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

%% Set-up
% intrinsic curvature constants
kc1 = 0.002;
kc2 = 0.004;

% arclength values
L = 90;
s_crit = 45;
ds = 0.5;
s = 0:ds:L;

%% Calculate the minimum kc value for 1-layer conversion
kc_optim = (kc2 * (3*L^2*s_crit - 3*L*s_crit^2 + (L - s_crit)^3 ) + kc1*s_crit^3) ...
         /(3*L^2*s_crit - 3*L*s_crit^2 + (L - s_crit)^3 + s_crit^3)

%% calculate 
k0_2 = kappa_0_2layer(s, kc1, kc2, L, s_crit);
k0_1_optim = kappa_0_1layer(s, kc_optim, L);

%% Plot the kappa_0 functions
fk0 = figure(1);
plot(s, k0_2, 'b-','DisplayName', 'actual 2-layer'); hold on;
plot(s, k0_1_optim, 'r-', 'DisplayName', 'optim 1-layer'); hold off;
xline(s_crit, 'r--', 'displayname', 'Tissue Boundary');
xlabel('s [mm]'); ylabel('\kappa_0 [1/mm]');
grid on;
legend();

%% Get the shapes
w_init_2layer = (kc1 * (s_crit/L)^2 + kc2 * (1 - s_crit/L)*(1 + s_crit/L)) * [1; 0 ;0];
[~, pmat_2layer, ~] = fn_intgEP_v3_2layers(w_init_2layer,kc1,kc2,s_crit,0,0,ds,length(s),B,Binv);

w_init_1layer_optim = [kc_optim; 0; 0];
[~, pmat_1layer_optim, ~] = fn_intgEP_v1_1layer(w_init_1layer_optim,kc_optim,0,0,ds,length(s),B,Binv);

%% Plot the 2-d shapes
f2d = figure(2);
plot(pmat_2layer(3,:), pmat_2layer(2,:), 'b-', 'displayname', 'actual 2-layer'); hold on;
plot(pmat_1layer_optim(3,:), pmat_1layer_optim(2,:), 'r-', 'displayname', 'optim 1-layer'); hold off;
xline(s_crit, 'r--', 'displayname', 'Tissue Boundary');
xlabel('z [mm]'); ylabel('y [mm]');
grid on;
legend()

%% Functions
% k0 for 1 layer case
function k0 = kappa_0_1layer(s, kc1, L)
    k0 = kc1 * (1 - s/L).^2;
    
end

% k0 for 2 layer case
function k0 = kappa_0_2layer(s, kc1, kc2, L, s_crit)
    % separate the arclength values
    s_before = s(s <= s_crit);
    s_after = s(s > s_crit);
    
    % calculate the curvature
    k0_before = kc1 * ( s_crit - s_before).^2 / L^2 ...
       + kc2 * (1 - s_crit/L)*(1 - 2*s_before/L + s_crit/L);
   
    k0_after = kc2 * (1 - s_after/L).^2;
    
    k0 = [k0_before k0_after];
    
end

