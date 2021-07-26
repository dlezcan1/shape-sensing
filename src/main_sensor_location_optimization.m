%% main_sensor_location_optimization
clear; close all;
%% Setup
% Needle shape generation params
params.kc = 0.0025508;
params.w_init = [];
params.sigma = 0.005;
params.L = 90;
params.ds = 0.5;
params.sloc_ds = params.ds*floor(2.5 / params.ds); % change only the numerator 
params.num_AA = 3;
params.divide_needle = false; % divide the needle equally in AAs?
params.fmincon = false; % whether to use fmincon or not

% cost functions for SLO
cost_fn_tip_mean     = @(slocs) cost_sensor_optimization(slocs, 'kc', params.kc, ...
                            'w_init', params.w_init, 'sigma', params.sigma, 'L', params.L, ...
                            'type', "tip-mean", 'ds', params.ds);
cost_fn_tip_bounds   = @(slocs) cost_sensor_optimization(slocs, 'kc', params.kc, ...
                            'w_init', params.w_init, 'sigma', params.sigma, 'L', params.L, ...
                            'type', "tip-bounds", 'ds', params.ds);
cost_fn_shape_bounds = @(slocs) cost_sensor_optimization(slocs, 'kc', params.kc, ...
                            'w_init', params.w_init, 'sigma', params.sigma, 'L', params.L, ...
                            'type', "shape-bounds", 'ds', params.ds);
cost_fn_shape_mean   = @(slocs) cost_sensor_optimization(slocs, 'kc', params.kc, ...
                            'w_init', params.w_init, 'sigma', params.sigma, 'L', params.L, ...
                            'type', "shape-mean", 'ds', params.ds);

% Saving name
save_bool = true;  % change if you'd like to save
data_dir = fullfile("../data/");
filesave_base = "SLO_numAA-%d";                        
               

%% Optimization
% initialiation
bnds = linspace(0, params.L, params.num_AA + 1)';
slocs_0 = mean([bnds(1:end-1), bnds(2:end)],2); % halfway points between sections
% slocs_0 = sort(randperm(params.L, params.num_AA))'; % random points along the needle

disp("Initial sensor locations:");
disp(slocs_0');

% bounds
if params.divide_needle
    lb = bnds(1:end-1);
    ub = bnds(2:end);
else
    lb = zeros(params.num_AA, 1);
    ub = params.L * ones(params.num_AA, 1);
end
lb(lb <= 0)        = lb(lb<= 0) + params.sloc_ds;
ub(ub >= params.L) = ub(ub >= params.L) - params.sloc_ds;
assert(all(lb < ub))

disp("bounds:");
disp([lb, ub]);

% fmincon optimization
if params.fmincon
    % options
    error("NotImplementedError: fmincon is not supported.");
    warning('off');
    disp("Beginning optimization...");
    tic; 
    maxiter = 1000;
    timelimit = 100;
    tolfun = 1e-6;%1e-4;

    oldopts = optimset('fmincon');
    %psoldopts = psoptimset;
    scalef = 1/cost_fn(slocs_0);
    Tol = 1e-14*10^(ceil(log10(scalef)));

    % options = optimset(oldopts,'Algorithm','interior-point','TolFun',Tol,'TolX',1e-8,...
    %     'MaxFunEvals',10000,'Display','notify','MaxIter',maxiter,'MaxFunEvals',100000);
    options = optimset(oldopts,'Algorithm','sqp','TolFun',Tol,'TolX',1e-2,...
        'Display','iter-detailed','MaxIter',maxiter,'MaxFunEvals',100000);

    % optimization
    [slocs_optim, fval, exitflag] = fmincon(@(sloc) cost_fn(sloc), ...
                                     slocs_0, [], [], [], [], lb, ub, [], options);

    toc
    warning('on');
    disp(" ");

% brute force optimization
else
    % gather all of the tests
    slocs_test = (lb(1):params.sloc_ds:ub(1))';
    for i = 2:params.num_AA
        slocs_test = cart_product(slocs_test, (lb(i):params.sloc_ds:ub(i))');
    end
    
    % remove unordered values
    mask_include = all(diff(slocs_test, 1, 2) > 0, 2);
    slocs_test = slocs_test(mask_include,:);
    
    % create data table
    data_tbl = array2table(slocs_test, 'VariableNames', "SensorLocation_" + (1:params.num_AA));
    
    % iterate through all of the costs
    warning('off');
    progressbar('Brute Force Optimization')
    for i = 1:size(data_tbl,1)
       data_tbl.cost_tip_mean(i)     = cost_fn_tip_mean(data_tbl{i,1:params.num_AA}');
       data_tbl.cost_tip_bounds(i)   = cost_fn_tip_bounds(data_tbl{i,1:params.num_AA}');
       data_tbl.cost_shape_mean(i)   = cost_fn_shape_mean(data_tbl{i,1:params.num_AA}');
       data_tbl.cost_shape_bounds(i) = cost_fn_shape_bounds(data_tbl{i,1:params.num_AA}');
       
       progressbar(i/size(data_tbl,1));
    end
    warning('on');
    
    % determine optimal slocs
%     [min_cost,min_cost_idx] = min(data_tbl.cost);
%     slocs_optim = data_tbl{min_cost_idx,1:params.num_AA};
end


% Display results
disp("Optimal sensor locations:")
disp(reshape(slocs_optim, 1, []));

%% Save results
if save_bool
    outfile = fullfile(data_dir, sprintf(filesave_base, params.num_AA) + ".mat");
    outfile = filename_no_overwrite(outfile);

    if params.fmincon
        save(outfile, 'params', 'lb', 'ub',...
            'slocs_optim', 'fval', 'exitflag', 'options');
    else
        save(outfile, 'params', 'lb', 'ub',...
            'data_tbl');
        writetable(data_tbl, strrep(outfile, '.mat', '.xlsx'));
        fprintf("Saved data table to: %s\n", strrep(outfile, '.mat', '.xlsx'));
    end
        
    fprintf("Saved data file: %s\n\n", outfile);
end