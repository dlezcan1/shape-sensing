%% main_sensor_location_optimization
clear; 
%% Setup
% Needle shape generation params
params.kc = 0.0025508;
params.w_init = [];
params.sigma = 0.005;
params.L = 90;
params.ds = 0.5;
params.num_AA = 3;
params.divide_needle = false; % divide the needle equally in AAs?

% cost function for SLO
weight_bounds = 1.0; % where to change
weight_mean = 1 - weight_bounds;
cost_fn = @(slocs) weight_bounds * cost_sensor_optimization(slocs, 'kc', params.kc, ...
                            'w_init', params.w_init, 'sigma', params.sigma, 'L', params.L, ...
                            'type', "tip-bounds", 'ds', params.ds) + ...
                   weight_mean   * cost_sensor_optimization(slocs, 'kc', params.kc, ...
                            'w_init', params.w_init, 'sigma', params.sigma, 'L', params.L, ...
                            'type', "tip-mean", 'ds', params.ds);
cost_fn2 = @(x) cost_fn(params.L*x);

% Saving name
save_bool = false;  % change if you'd like to save
data_dir = fullfile("../data/");
filesave_base = fullfile(data_dir, "SLO_numAA-%d");                        
               

%% Optimization
% initialiation
bnds = linspace(0, 1, params.num_AA + 1)';
slocs_0 = mean([bnds(1:end-1), bnds(2:end)],2); % halfway points between sections
slocs_0 = sort(rand(params.num_AA,1));
disp("Initial sensor locations:")
disp(params.L * slocs_0')


% bounds
if params.divide_needle
    
    lb = bnds(1:end-1);
    ub = bnds(2:end);
else
    lb = zeros(params.num_AA, 1);
    ub = 1 * ones(params.num_AA, 1);
end

disp("bounds:");
disp([lb, ub]);

% options
warning('off');
disp("Beginning optimization...");
tic; 
maxiter = 1000;
timelimit = 100;
tolfun = 1e-6;%1e-4;

oldopts = optimset('fmincon');
%psoldopts = psoptimset;
scalef = 1/cost_fn2(slocs_0);
Tol = 1e-14*10^(ceil(log10(scalef)));

% options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
%     'MaxFunEvals',10000,'Display','notify','MaxIter',maxiter,'MaxFunEvals',100000);
options = optimset(oldopts,'Algorithm','sqp','TolFun',Tol,'TolX',1e-2,...
    'Display','iter-detailed','MaxIter',maxiter,'MaxFunEvals',100000);

% optimization
[slocs_optim, fval, exitflag] = fmincon(@(sloc) cost_fn2(sloc), ...
                                 slocs_0, [], [], [], [], lb, ub, [], options);

toc
warning('on');
disp(" ");

% Display results
disp("Optimal sensor lcoations:")
disp(reshape(params.L*slocs_optim, 1, []));

%% Save results
if save_bool
    outfile = sprintf(filesave_base, params.num_AA) + ".mat";
    outfile = filename_no_overwrite(outfile);

    save(outfile, 'params', 'lb', 'ub', 'weight_bounds', 'weight_mean'); 
    fprintf("Saved data file: %s\n\n", outfile);
end