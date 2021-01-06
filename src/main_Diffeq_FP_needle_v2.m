%% main_Diffeq_FP_needle.m
%
% Solve the F-P equation for the prob. density of deformation
%
% - written by Dimitri Lezcano

global kc L w dw s_l w_i w_j w_k ds sigma A0 G0

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

% insertion parameters
L = 90;
kc = 0.005;
sigma = 0.0001; % gaussian noise uncertainty

% arclength coordinate
ds = 1;
s_l = 0:ds:L;
N = length(s_l);

% angular deformation coordinates
dw = 0.0005;
w = -.01:dw:.02;
M = length(w);

w_i = w; % omega_x vector
w_j = w; % omega_y vector
w_k = w; % omega_z vector


%% boundary conditions

% initial deformation distribution
w_init = [kc; 0; 0];
[~, idx_i_init] = min(abs(w - w_init(1)));
[~, idx_j_init] = min(abs(w - w_init(2)));
[~, idx_k_init] = min(abs(w - w_init(3)));

p_init = zeros(M, M, M);
p_init(idx_i_init, idx_j_init, idx_k_init) = 1;


%% main loop
% solution tensor
prob = zeros(M, M, M, N); % w_x, w_y, w_z, s | i, j, k, l
prob = init_probability();
h = waitbar(0, 'Please wait...');
for l = 2:N               % arclength coordinate
    for i = 2:M-1         % w1 | omega_x
        for j = 2:M-1     % w2 | omega_y
            for k = 2:M-1 % w3 | omega_z
                
                % system matrix
                K_ijk = gen_upmatrix(i, j, k, l);
                
                % local probability vector
                p_vect = [prob(i-1, j, k, l); prob(i+1, j, k, l);
                          prob(i, j-1, k, l); prob(i, j+1, k, l);
                          prob(i, j, k-1, l); prob(i, j, k+1, l);];
                
                % update the point      
                prob(i,j,k,l) = -1/ds* prob(i,j,k,l-1) + K_ijk * p_vect;
                prob(i,j,k,l) = -prob(i,j,k,l)/(1/ds + 3*sigma^2/dw^2);
                
            end
        end
    end   
    
    waitbar(l/N, h, sprintf('Simulation: %.2f%% Complete.', l/N*100));
   
end
waitbar(l/N, h, 'Performing normalization.');
prob = prob./sum(prob, [1,2,3]); % normalize the probability densities
close(h);

%% plots
f1 = figure(1);
set(gcf,'units','normalized','position', [1/6, 1/6, 4/6, 4/6]);

for l = 1:N
    prob_l = prob(:,:,:,l);
    [w1_grid, w3_grid] = meshgrid(w_i, w_k);
    ind = find(max(prob_l,[],'all') == prob_l);
    [idx_i, idx_j, idx_k] = ind2sub(size(prob_l), ind);
    w1_max = w_i(idx_i); w2_max = w_j(idx_j); w3_max = w_j(idx_k);
    
    subplot(1,2,1);
    surf(w3_grid, w1_grid, reshape(prob_l(:,idx_j,:),[M M]));
    xlabel('\omega_z'); ylabel('\omega_x'); zlabel('probability');
    title('\omega_z & \omega_x') 
    zlim([0 1.1*max(prob_l,[],'all')])
    view([ -57.0300    9.8047])
    
    subplot(1,2,2);
    [w2_grid, w3_grid] = meshgrid(w_j, w_k);
    surf(w3_grid, w2_grid, reshape(prob_l(idx_i,:,:),[M M]));
    xlabel('\omega_z'); ylabel('\omega_y'); 
    title('\omega_z & \omega_y') 
    zlim([0 1.1*max(prob_l,[],'all')])
    view([ -57.0300    9.8047])
    
    [w1_grid, w2_grid] = meshgrid(w_i, w_j);
    surf(w3_grid, w2_grid, reshape(prob_l(:,:,idx_k),[M M]));
    xlabel('\omega_x'); ylabel('\omega_y'); 
    title('\omega_z & \omega_y') 
    zlim([0 1.1*max(prob_l,[],'all')])
    view([ -57.0300    9.8047])
    
    
    sgtitle(sprintf("s = %.4f mm | %s_{max} = [%f, %f, %f]", ...
        s_l(l), "\omega",w1_max, w2_max, w3_max));
    pause(.25);
    
end

%% Functions

% kappa_0 functions
function k0 = kappa_0(s)
    global kc L
    
    k0 = kc * (1 - s/L).^2;
    
end

function k0_prime = kappa_0_prime(s)
    global kc L
    
    k0_prime = -2 * kc * (1 - s/L);

end

% diff-eq. update functions
function [ K_ijk, K_i, K_j, K_k ] = gen_upmatrix(i, j, k, l)
% Function to generate the matrices for the update
% Input
%   indices i,j,k,l
%
% Return
%   K_ijk, the update row matrix combined [K_i, K_j, K_k]
%   K_i, K_j, K_k, the update matrices for index (i, j, k, l)


    global ds dw sigma A0 G0 s_l w_i w_j w_k
    
    % parameter retrieval
    s = s_l(l);
    w1 = w_i(i);
    w2 = w_j(j);
    w3 = w_k(k);
    
    % kappa_0 retrieval
    k0 = kappa_0(s);
    k0_p = kappa_0_prime(s);
    
    % the update matrices for i, j, k
    K_i = [ -(k0_p + (A0-G0)/A0 * w2 .* w3)/(2*dw) - sigma^2/dw^2, ...
             (k0_p + (A0-G0)/A0 * w2 .* w3)/(2*dw) - sigma^2/dw^2];
    
    K_j = [ -(k0 .* w3 + (G0-A0)/A0 .* w1 .* w3)/(2*dw) - sigma^2/dw^2 ...
             (k0 .* w3 + (G0-A0)/A0 .* w1 .* w3)/(2*dw) - sigma^2/dw^2];
    
    K_k = [ A0/G0 * k0 .* w2/(2*dw) - sigma^2/dw^2 ...
            -A0/G0 * k0 .* w2/(2*dw) - sigma^2/dw^2 ];
        
    K_ijk = [ K_i, K_j, K_k ];
        
    
    % division term (if needed, can be implemented in the loop)
    div = -(1/ds + 3*sigma^2/dw^2);
    
end

function p_init = init_probability()
% generate dirac delta from omega_0(s) first.
    global s_l w
    
    k0 = kappa_0(s_l);
    
    M = length(w);
    N = length(s_l);
    
    p_init = zeros(M, M, M, N);
    
    for l = 1:1
        w_0 = [k0(l); 0; 0];
        
        [~, idx_i_init] = min(abs(w - w_0(1)));
        [~, idx_j_init] = min(abs(w - w_0(2)));
        [~, idx_k_init] = min(abs(w - w_0(3)));

        p_init(idx_i_init, idx_j_init, idx_k_init, l) = 1;
        
    end
    
end

    
    