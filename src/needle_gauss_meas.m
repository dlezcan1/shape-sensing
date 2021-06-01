function [pos, pos_bounds, wv, Sigma_wv] = needle_gauss_meas(slocs, parameters)
    arguments
        slocs (1,:);
        parameters.kc double = 0.0025508;
        parameters.w_init = [];
        parameters.sigma double = 0.0025;
        parameters.L double = 90;
    end
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
    L = parameters.L;
    kc = parameters.kc; %0.003;
    sigma = parameters.sigma; % gaussian noise uncertainty
    w_init = parameters.w_init;

    % arclength coordinate
    ds = 0.5;

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
    S.kc = parameters.kc;
    S.sigma = sigma; % standard deviation
    S.sigma2 = S.sigma^2;
    if isempty(w_init)
        S.w_init = [kc; 0; 0];
    else
        S.w_init = w_init;
    end
    S.ds = ds;
    S.arclengths = 0:S.ds:S.L;

    %- measurement parameters
    S.cov_m = [0.0520, 0.0532, inf]; % diagonal covariance of measurement probability
    S.meas_locations = slocs;
    [~, S.meas_idxs] = min(abs(S.arclengths' - S.meas_locations), [], 1); % indexes just in case
    S.appx_meas_locations = S.arclengths(S.meas_idxs); % approximate locations (if not already there)
    
    %% Break up arclengths for each measurement cycle
    S.arclengths_per_measurement = {};
    idx_0 = 1;
    for i = 1:length(S.meas_locations)
       idxs = idx_0:S.meas_idxs(i);
       S.arclengths_per_measurement{i} = S.arclengths(idxs);
       idx_0 = S.meas_idxs(i);
    end
    S.arclengths_per_measurement{i+1} = S.arclengths(idx_0:end);
    
    %% Dynamical formulation
    needle_f = @(s, x) gaussian_dynamics(s, x, S);
    x0 = [S.w_init; 0e-8*ones(3,1)];

    % iterate through each of the measurements
    x_init = x0;
    s = []; x = [];
    for i = 1:length(S.arclengths_per_measurement)
       % iterate up to next AA location
       [s_i, x_i] = ode45(@(s,x) needle_f(s,x), S.arclengths_per_measurement{i}, x_init);

       % Bayesian fusion from ideal sensor measurement
       if ismember(s_i(end), S.appx_meas_locations)
            xp_AAi = posterior_measurement_ideal(s_i(end), x_i(end,:), S);
            x_i(end,:) = reshape(xp_AAi, 1, []);
       end

       % append new data
       if i < length(S.arclengths_per_measurement)
           s = [s; s_i(1:end-1)];
           x = [x; x_i(1:end-1,:)];
       else
           s = [s; s_i];
           x = [x; x_i];
       end
       % update x_init
       x_init = reshape(x_i(end,:), [], 1);
    end

    wv = x(:,1:3)';
    sigma2_vect = x(:,4:6)';
    Sigma_wv = zeros(3,3,size(wv, 2));
    for i = 1:size(Sigma_wv, 3)
       Sigma_wv(:,:,i) = diag(sigma2_vect(:, i)); 
    end
    
    %% Generate the probabilistic needle shapes
    % mean shape
    w0 = wv;
    w0_prime = [diff(w0, 1, 2), zeros(3,1)];
    [~,pos,~] = fn_intgEP_w0_Dimitri(S.w_init, w0, w0_prime, 0, 0, ds, size(wv, 2), B, Binv);

    % sample Gaussian omega (w3 is assumed to be fixed)
    theta_vals = linspace(0, 2*pi-0.001, 50);
    u = [cos(theta_vals); sin(theta_vals)];
    w_err_mat = zeros(length(theta_vals), 3, size(Sigma_wv, 3));

    for i = 1:length(theta_vals)
        u_i = u(:,i);
        for l = 1:size(Sigma_wv, 3)
            sig_i = Sigma_wv(1:2,1:2,l); % only grab first two

            w12_i = sig_i * u_i + wv(1:2,l);

            w_err_mat(i, 1:3, l) = [w12_i; wv(3,l)];

        end
    end 

    % shape bounds
    pos_bounds = zeros(size(w_err_mat));
    for i = 1:size(w_err_mat, 1)
        % grab w0 values
        w0_i = squeeze(w_err_mat(i,:,:));

        % approximate w0_prime
        w0_prime_i = [diff(w0_i, 1, 2), zeros(3,1)];

        % generate sampled shape
        [~,pmat,~] = fn_intgEP_w0_Dimitri(S.w_init, w0_i, w0_prime_i, 0, 0, ds, size(pos_bounds, 3), B, Binv);

        pos_bounds(i,:,:) = pmat;

    end
    
end

%% Helper Functions
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

% calculate posterior measurement update based on ideal deformation
function [xp, mu_p, cov_p] = posterior_measurement_ideal(s, x, S)
    % unpack x
    mu = reshape(x(1:3), 3,1);
    cov = reshape(x(4:6), 3,1);

    % generate measurement mean
    mu_m = kappa_0(s, S) * [1; 0; 0];
    cov_m = S.cov_m;

    [mu_p, cov_p] = gauss_bayes_diagonal(mu, cov, mu_m, cov_m);
    xp = [mu_p; cov_p];

end