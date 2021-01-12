%% main_Diffeq_FP_needle.m
%
% Solve the F-P equation for the prob. density of deformation
%
% - written by Dimitri Lezcano

%% Variable set-up
clear; clc;

global kc L w dw s_l w_i w_j w_k ds sigma A0 G0 deletion_indices

%% saving options 
save_bool = true;

directory = "Data/";
file_base = directory + "DiffEq_Results_sigma_%.4f";

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
% w_init = [ 0.0035703; 0.00072161; -0.0086653 ]; % data

if isempty(w_init)
    file_base = file_base + "_ideal";
else
    file_base = file_base + "_data";
end

% arclength coordinate
ds = 1;
s_l = 0:ds:L;
N = length(s_l)

% angular deformation coordinates
dw = 0.002;
w = -.05:dw:.05;
M = length(w)

w_i = w; % omega_x vector
w_j = w; % omega_y vector
w_k = w; % omega_z vector

M_i = length(w_i); M_j = length(w_j); M_k = length(w_k); % length of each of these vectors

%% Indices for subset (removing boundary)
% find columns and rows to delete (indices)
sz = [M_i M_j M_k];
idxs_i = zeros(1, 2*M_j*M_k, 'uint32'); % i = 1 | M_i
for j = 1:M_j
    for k = 1:M_k
        idx_idx = sub2ind([M_j, M_k], j, k);
        idxs_i(idx_idx) = sub2ind(sz, 1, j, k); % i = 1
        idxs_i(idx_idx + M_j*M_k) = sub2ind(sz, M_i, j, k); % i = M_i

    end
end
clear i j k;

idxs_j = zeros(1, 2*M_i*M_k, 'uint32'); % j = 1 | M_j
for i = 1:M_i
    for k = 1:M_k
        idx_idx = sub2ind([M_i, M_k], i, k);
        idxs_j(idx_idx) = sub2ind(sz, i, 1, k); % j = 1
        idxs_j(idx_idx + M_i*M_k) = sub2ind(sz, i, M_j, k); % j = M_j

    end
end
clear i j k;

idxs_k = zeros(1, 2*M_i*M_j, 'uint32'); % k = 1 | M_k
for i = 1:M_i
    for j = 1:M_j
        idx_idx = sub2ind([M_i, M_j], i, j);
        idxs_k(idx_idx) = sub2ind(sz, i, j, 1); % k = 1
        idxs_k(idx_idx + M_i*M_j) = sub2ind(sz, i, j, M_k); % k = M_k

    end
end
clear i j k;

deletion_indices = unique([idxs_i idxs_j idxs_k]);

%% boundary conditions
% initial boundary condition
prob = init_probability(w_init); % w1, w2, w3, s | wx, wy, wz, s | i, j, k, l

%% Test gen_sysmatrix2
% disp('Gen 1');
% tic;
% K1 = gen_sysmatrix(1);
% toc
% 
% disp('Gen 1 Parallel');
% parpool('local')
% tic;
% K2 = gen_sysmatrix(1);
% toc
% 
% dK = full(abs(K1 - K2));
% fprintf('Max: %.5f, Mean: %.5f\n', max(dK, [], 'all'), mean(dK, 'all'));
% return;

%% main loop
h = waitbar(0, 'Please wait...', 'units', 'normalized', 'position', [.4 .455 .2 .09]);
cum_time = 0; % cumulative time for averaging
for l = 2:N               % arclength coordinate
    tic;
    % system matrix
    [~, K_sub] = gen_sysmatrix(l);
%     waitbar((l-1)/N, h, sprintf('Generated K matrix: %.2f%% Complete.', l/N*100))

    % vectorize 
    prob_lm1_vect = prob_vectorize(prob(2:M-1,2:M-1,2:M-1,l-1), 1); % p_{2:M-1,2:M-1,2:M-1,l-1}

    % update new vector
    prob_l_vect = K_sub\prob_lm1_vect; % regular least squares
%     prob_l_vect = lsqnonneg(K_sub, prob_lm1_vect); % non-negative least-squares
    
    % update the new array
    prob(2:M-1,2:M-1,2:M-1,l) = prob_tensorize(prob_l_vect, M-2, 1);
    
    % time calculations
    cum_time = cum_time + toc;
    mean_time = cum_time/(l-1); % per run
    approx_time_left = (N-l)*mean_time;
    
    waitbar(l/N, h, sprintf('Simulation: %.2f%% Complete.\n Time per iteration: %.2fs \nEstimated time to complete: %d min %d secs',...
        l/N*100, mean_time, floor(approx_time_left/60), ceil(mod(approx_time_left,60) )));

end

waitbar(l/N, h, 'Performing normalization.');
prob = prob./sum(prob, [1,2,3]); % normalize the probability densities
close(h);

%% save the current run
if save_bool
    file_base = sprintf(file_base, sigma);
    save(file_base + ".mat");
    disp("Saved: " + file_base + ".mat")
end

%% plots
f1 = figure(1);
set(gcf,'units','normalized','position', [1/20, 1/20, 18/20, 18/20]);
w_interp = linspace(min(w), max(w), 150); % for interpolation
if save_bool
    vidfile = VideoWriter(file_base + "_results.mp4",'MPEG-4');
    vidfile.FrameRate = 2;
    open(vidfile);
end

std_filt = 4;

for l = 1:N
    if l == 1
        z_lim = 0.2;
        
    elseif l > 5
        z_lim = 0.02;
        
    end
    
    % find maximum indices
    prob_l = prob(:,:,:,l);
    z_lim = max(prob_l,[],'all');

    ind = find(max(prob_l,[],'all') == prob_l);
    [idx_i, idx_j, idx_k] = ind2sub(size(prob_l), ind);
    w1_max = w_i(idx_i); w2_max = w_j(idx_j); w3_max = w_k(idx_k);
        
    % x-y plots
    subplot(2,3,1);    
    [w1_grid, w2_grid] = meshgrid(w_i, w_j);
    prob_grid = reshape(prob_l(:,:,idx_k), [M_i, M_j]);

    surf(w1_grid, w2_grid, prob_grid, 'edgecolor', 'none');
%     contour3(w1_grid, w2_grid, prob_grid, 'linewidth', 2);
    xlabel('\omega_y'); ylabel('\omega_x'); 
    title('\omega_x & \omega_y') 
    zlim([0 z_lim])
    
    subplot(2,3,4);
    surf(w1_grid, w2_grid, prob_grid, 'edgecolor', 'none');
%     contour3(w1_grid, w2_grid, prob_grid, 'linewidth', 2);
    xlabel('\omega_y'); ylabel('\omega_x'); 
    title('\omega_x & \omega_y') 
    zlim([0 z_lim])
    view([ 0    90])
       
    % x-z plots
    subplot(2,3,2);
    [w1_grid, w3_grid] = meshgrid(w_i, w_k);
    prob_grid = reshape(prob_l(:,idx_j,:), [M_i, M_k]);

    surf(w1_grid, w3_grid, prob_grid, 'edgecolor', 'none');
%     contour3(w1_grid, w3_grid, prob_grid, 'linewidth', 2);
    xlabel('\omega_z'); ylabel('\omega_x'); zlabel('probability');
    title('\omega_z & \omega_x') 
    zlim([0 z_lim])
    
    subplot(2,3,5);
    surf(w1_grid, w3_grid, prob_grid, 'edgecolor', 'none');
%     contour3(w1_grid, w3_grid, prob_grid, 'linewidth', 2);
    xlabel('\omega_z'); ylabel('\omega_x'); zlabel('probability');
    title('\omega_z & \omega_x') 
    zlim([0 z_lim])
    view([ 0    90])
    
    % y-z plots
    subplot(2,3,3);
    [w2_grid, w3_grid] = meshgrid(w_j, w_k);
    prob_grid = reshape(prob_l(idx_i,:,:), [M_j, M_k]);

    surf(w2_grid, w3_grid, prob_grid, 'edgecolor', 'none');
%     contour3(w2_grid, w3_grid, prob_grid, 'linewidth', 2);
    xlabel('\omega_z'); ylabel('\omega_y'); 
    title('\omega_z & \omega_y') 
    zlim([0 z_lim])
    
    
    subplot(2,3,6);
    surf(w2_grid, w3_grid, prob_grid, 'edgecolor', 'none');
%     contour3(w2_grid, w3_grid, prob_grid, 'linewidth', 2);
    xlabel('\omega_z'); ylabel('\omega_y'); 
    title('\omega_z & \omega_y') 
    zlim([0 z_lim])
    view([ 0    90])    
    
    sgtitle(sprintf("s = %.4f mm | %s_{max} = [%f, %f, %f]", ...
        s_l(l), "\omega",w1_max, w2_max, w3_max));
    if save_bool
        F(ind) = getframe(gcf); 
        writeVideo(vidfile,F(ind));
        
    end
    
    pause(.05);
    
end

if save_bool
    close(vidfile);
    disp("Wrote video file: " + file_base + "_results.mp4");
end

%% Functions

% kappa_0 functions
function k0 = kappa_0(s)
    global kc L
    
    k0 = kc * (1 - s/L).^2;
    
end

function k0_prime = kappa_0_prime(s)
    global kc L
    
    k0_prime = -2 * kc/L * (1 - s/L);

end

function p_init = init_probability(w_init)
% generate dirac delta from omega_0(0) first.
    global s_l w
    
    if isempty(w_init)
        k0 = kappa_0(s_l);
        w_init = [k0(1); 0; 0];
    end

    M = length(w);
    N = length(s_l);

    p_init = zeros(M, M, M, N);

    for l = 1:1
        [~, idx_i_init] = min(abs(w - w_init(1)));
        [~, idx_j_init] = min(abs(w - w_init(2)));
        [~, idx_k_init] = min(abs(w - w_init(3)));

        p_init(idx_i_init, idx_j_init, idx_k_init, l) = 1;

    end
    
end

function [K, K_sub] = gen_sysmatrix(l)
% generate the system matrix K * u = f
%
% Args:
%   prob is the M x M x M x N tensor
%   l is the arclength index
    global w_i w_j w_k s_l ds dw sigma A0 G0 deletion_indices
    % Iniital setup
    M_i = length(w_i); M_j = length(w_j); M_k = length(w_k); % variable just in case
    sz = [M_i, M_j, M_k]; % size of prob. matrix
    
    k0 = kappa_0(s_l(l));
    k0_p = kappa_0_prime(s_l(l));
    
    K = sparse(M_i * M_j * M_k, M_i * M_j * M_k);
    
    for k = 2:M_k-1
        for j = 2:M_j-1
            for i = 2:M_i-1
                % flattened local indices
%                 indx_ijk = M_j * M_i * (k-1) + M_i * (j-1) + (i-1) + 1; %i,j,k
%                 indx_im = M_j * M_i * (k-1) + M_i * (j-1) + (i-2) + 1; % i-1 
%                 indx_ip =  M_j * M_i * (k-1) + M_i * (j-1) + (i) + 1; % i+1
%                 indx_jm =  M_j * M_i * (k-1) + M_i * (j-2) + (i-1) + 1;   % j-1
%                 indx_jp =  M_j * M_i * (k-1) + M_i * (j) + (i-1) + 1; % j+1
%                 indx_km =  M_j * M_i * (k-2) + M_i * (j-1) + (i-1) + 1; % k-1
%                 indx_kp =  M_j * M_i * (k) + M_i * (j-1) + (i-1) + 1; % k+1
                
                indx_ijk = sub2ind([M_i, M_j, M_k], i, j, k); % i,j,k
                indx_im = sub2ind([M_i, M_j, M_k], i-1, j, k); % i-1
                indx_ip = sub2ind([M_i, M_j, M_k], i+1, j, k); % i+1
                indx_jm = sub2ind([M_i, M_j, M_k], i, j-1, k); % j-1
                indx_jp = sub2ind([M_i, M_j, M_k], i, j+1, k); % j+1
                indx_km = sub2ind([M_i, M_j, M_k], i, j, k-1); % k-1
                indx_kp = sub2ind([M_i, M_j, M_k], i, j, k+1); % k+1
                
                % get omega values
                w1 = w_i(i); w2 = w_j(j); w3 = w_k(k);
                
                % update row_ijk with the corresponding values
                K(indx_ijk, indx_ijk) = 1/ds + 3 * sigma^2/dw^2;
                
                % (i-1) & (i+1)
                K(indx_ijk, indx_im) = 1/2 * (-(k0_p + (A0-G0)/A0 * w2 * w3)/(dw) ...
                    - sigma^2/dw^2);
                K(indx_ijk, indx_ip) = 1/2 * ((k0_p + (A0-G0)/A0 * w2 * w3)/(dw) ...
                    - sigma^2/dw^2);
                
                % (j-1) & (j+1)
                K(indx_ijk, indx_jm) = 1/2 * (-(k0 * w3 + (G0-A0)/A0 * w1 * w3)/(dw) ...
                    - sigma^2/dw^2);
                K(indx_ijk, indx_jp) = 1/2 * ((k0 * w3 + (G0-A0)/A0 * w1 * w3)/(dw) ...
                    - sigma^2/dw^2);
                
                % (k-1) & (k+1)
                K(indx_ijk, indx_km) = 1/2 * (A0/G0 * k0 * w2/(dw) - sigma^2/dw^2);
                K(indx_ijk, indx_kp) = 1/2 * (-A0/G0 * k0 * w2/(dw) - sigma^2/dw^2);
                
            end
        end
    end
    
    K = ds * K; % mutiply the division term
    K_sub = K;
    
    % delete the rows and columns for the boundary points
    K_sub(deletion_indices,:) = []; % rows to delete
    K_sub(:,deletion_indices) = []; % cols to delete
    
end

function [K, K_sub] = gen_sysmatrix2(l)
% generate the system matrix K * u = f
%
% vectorized version
%
% Args:
%   prob is the M x M x M x N tensor
%   l is the arclength index
    global w_i w_j w_k s_l ds dw sigma A0 G0 deletion_indices
    % Iniital setup
    M_i = length(w_i); M_j = length(w_j); M_k = length(w_k); % variable just in case
    sz = [M_i, M_j, M_k]; % size of prob. matrix
    
    k0 = kappa_0(s_l(l));
    k0_p = kappa_0_prime(s_l(l));
    
    K = sparse(M_i * M_j * M_k, M_i * M_j * M_k);
    
    % index values
    I = (2:M_i-1)';
    J = (2:M_j-1)';
    K_idx = (2:M_k-1)';   
    indxs = cart_product(I, J, K_idx); % matrix of index positions

    
    % prepare omega vectors values
    w1 = reshape(w_i, [], 1);
    w2 = reshape(w_j, [], 1);
    w3 = reshape(w_k, [], 1);
    
    W1 = w1(indxs(:,1)); % vectorized omega1 alighed with indxs
    W2 = w2(indxs(:,2)); % vectorized omega2 aligned with indxs
    W3 = w3(indxs(:,3)); % vectorized omega3 aligned with indxs
    
    % diagonal elements
    indxs_ijk = sub2ind([M_i, M_j, M_k], indxs(:,1), indxs(:,2), indxs(:,3));
    K(sub2ind(size(K), indxs_ijk, indxs_ijk)) = ...
        (1/ds + 3*sigma^2/dw^2);
    
    
    % index: i
    indxs_tmp = sub2ind([M_i, M_j, M_k], indxs(:,1)-1, indxs(:,2), indxs(:,3)); % i-1
    K(sub2ind(size(K), indxs_ijk, indxs_tmp)) = ... 
        1/2 * (-(k0_p + (A0-G0)/A0 * W2 .* W3)/(dw) - sigma^2/dw^2); % i-1
    
    indxs_tmp = sub2ind([M_i, M_j, M_k], indxs(:,1)+1, indxs(:,2), indxs(:,3)); % i+1
    K(sub2ind(size(K), indxs_ijk, indxs_tmp)) = ...
        1/2 * ((k0_p + (A0-G0)/A0 * W2 .* W3)/(dw) - sigma^2/dw^2); % i+1
    
    % index: j
    indxs_tmp = sub2ind([M_i, M_j, M_k], indxs(:,1), indxs(:,2)-1, indxs(:,3)); % j -1
    K(sub2ind(size(K), indxs_ijk, indxs_tmp)) = ... 
        1/2 * (-(k0 * W3 + (G0-A0)/A0 * W1 .* W3)/(dw) - sigma^2/dw^2); % j-1
    
    indxs_tmp = sub2ind([M_i, M_j, M_k], indxs(:,1), indxs(:,2)+1, indxs(:,3)); % j + 1
    K(sub2ind(size(K), indxs_ijk, indxs_tmp)) = ... 
        1/2 * ((k0 * W3 + (G0-A0)/A0 * W1 .* W3)/(dw) - sigma^2/dw^2); % j+1
    
    % index: k
    indxs_tmp = sub2ind([M_i, M_j, M_k], indxs(:,1), indxs(:,2), indxs(:,3)-1); % k - 1
    K(sub2ind(size(K), indxs_ijk, indxs_tmp)) = ... 
        1/2 * (A0/G0 * k0 * W2/(dw) - sigma^2/dw^2); % k-1
    
    indxs_tmp = sub2ind([M_i, M_j, M_k], indxs(:,1), indxs(:,2), indxs(:,3)+1); % k + 1
    K(sub2ind(size(K), indxs_ijk, indxs_tmp)) = ... 
        1/2 * (-A0/G0 * k0 * W2/(dw) - sigma^2/dw^2); % k+1
    
    clear indxs_tmp indxs_ijk
    
    % scale by ds
    K = ds * K; % multiply the division term
    K_sub = K;
    
    % delete the rows and columns for the boundary points
    K_sub(deletion_indices,:) = []; % rows to delete
    K_sub(:,deletion_indices) = []; % cols to delete
    
end

function p_vect = prob_vectorize(prob_tensor, N)
% function to vectorize probability tensor
% i -> j -> k
% for l = 1:
%   for k = 1:
%       for j = 1:
%           for i = 1;
%               vector = tensor(i,j,k,l);

    p_vect = reshape(prob_tensor, [], N);
    
end

function p_tens = prob_tensorize(prob_vector, M, N)
% function to vectorize probability tensor

    p_tens = reshape(prob_vector, M, M, M, N);
    
end

function [x_out, x1_new, x2_new, x3_new] = gaussian_smoothing(x, x1, x2, x3, sigma)
    % works only for 2-D
    
    % double the space for higer resolution
    x1_new = linspace(min(x1),max(x1),2*length(x1)-1);
    x2_new = linspace(min(x2),max(x2),2*length(x2)-1);
    x3_new = linspace(min(x3),max(x3),2*length(x3)-1);
    [x1_new_grid, x2_new_grid, x3_new_grid] = meshgrid(x1_new, x2_new, x3_new); 
    
    % start adding the peak values as gaussians
    x_out = zeros(2*size(x)-1); % 2-D
    for k = 1:length(x3)
        for j = 1:length(x2)
            for i = 1:length(x1)
                x1i = x1(i); x2j = x2(j); x3k = x3(k);
                A = x(i, j, k); % the i,j,k-th value of the peak
                dx2_grid = (x1_new_grid - x1i).^2 + (x2_new_grid - x2j).^2 + (x3_new_grid - x3k).^2;

                x_out = x_out + A*exp(-dx2_grid./(2*sigma^2));

            end
        end
    end
    x_out = x_out/sum(x_out, 'all');
end
            
            
    
    
