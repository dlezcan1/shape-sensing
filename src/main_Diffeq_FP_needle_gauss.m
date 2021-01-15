%% main_Diffeq_FP_needle_gauss.m
%
% Solve the F-P equation for the prob. density of deformation using Gaussian
% assumption
%
% - written by Dimitri Lezcano

%% Variable set-up
clear; clc;

%% saving options 
save_bool = false;

directory = "Data/";
file_base = directory + "DiffEq_Results_sigma_%.4f_gauss";

%% preamble
% physical parameters
Emod = 200e9*1e-6; % 200 GPa, conversion from N/m^2 to N/mm^2
Pratio = 0.29; % Poisson's ratio
diam = 0.9; % in mm
Ibend = pi*diam^4/64;

Gmod = Emod/2/(1+Pratio);
Jtor = pi*diam^4/32;

A0 = Emod*Ibend; 
G0 = Gmod*Jtor; 

B = diag([A0, A0, G0]);
Binv = inv(B);

% insertion parameters
L = 90;
kc = 0.0025508; %0.003;
sigma = 2.5*0.001; % gaussian noise uncertainty
w_init = []; % ideal case

if isempty(w_init)
    file_base = file_base + "_ideal";
else
    file_base = file_base + "_data";
end

% arclength coordinate
ds = 1;
s_l = 0:ds:L;
N = length(s_l);

% system parameters
%- constants
S.A0 = A0;
S.G0 = G0;
S.alpha1 = (S.A0 - S.G0)/S.A0;
S.alpha2 = - S.alpha1;
S.alpha3 = -S.A0 / S.G0;
S.B = B;
S.Binv = Binv;

%- insertion parameters
S.L = L;
S.kc = kc;
S.sigma = sigma; % standard deviation
S.sigma2 = S.sigma^2;
if isempty(w_init)
    S.w_init = [kc; 0; 0];
else
    S.w_init = w_init;
end
S.ds = ds;

%% Dynamical formulation
needle_f = @(s, x) gaussian_dynamics(s, x, S);
x0 = [S.w_init; 0e-8*ones(3,1)];
[s, x] = ode45(@(s, x) needle_f(s, x), [0, S.L], x0);

mu = x(:,1:3)';
sigma2_vect = x(:,4:6)';

%% Generate the probabilistic needle shapes

%% save the current run
if save_bool
    file_base = sprintf(file_base, sigma);
    save(file_base + ".mat");
    disp("Saved: " + file_base + ".mat")
end

%% plots

%% Functions
% kappa_0 functions
function k0 = kappa_0(s, S)
    k0 = S.kc * (1 - s/S.L).^2;
    
end

function k0_prime = kappa_0_prime(s, S)
    k0_prime = -2 * S.kc/S.L * (1 - s/S.L);

end

% mean and covariance dynamics
function dx = gaussian_dynamics(s, x, S)
% Args:
%   s = the arclength currently at
%   x = (mu; sigma^2) where mu and sigma^2 are the current dynamics

    % unpack x
    mu = x(1:3);
    sigma2 = x(4:6);
    
    % prepare dynamics
    k0 = kappa_0(s, S);
    k0_prime = kappa_0_prime(s, S);
    
    % dynamics formulas
    dmu = [k0_prime + S.alpha1 * mu(2) * mu(3);
           k0 * mu(3) + S.alpha2 * mu(1) * mu(3);
           S.alpha3 * k0 * mu(2)];
       
    dsigma2 = S.sigma2 * ones(3,1);
    
    % concatenate dynamics
    dx = [dmu; dsigma2];

end
