%% main_Diffeq.m
%
% Solve the Diffusion equation which describes the diffusion along the
% surronding tissue of sweat duct.
%
% - written by Jin Seob Kim

clear all

% this is useful when you want font size in your plots bigger
% set(0,'DefaultAxesFontSize',28); 

%% preamble
% physical parameters
D = 1; % Diffusion coefficient

% temporal coordinate
t_end = 10;
M = 101;
dt = t_end/(M-1);
t = [0:dt:t_end];

% radial coordinate
r_o = 10; % outer radius
r_i = 3; % inner radius

N = 101;
dr = (r_o - r_i)/(N-1);

r = [r_i:dr:r_o]';

%% boundary conditions
% boundary concentration (r = r_o)
u_out = 0;
u_in_0 = 10;

% boundary condition (r = r_i)
u_in = u_in_0*exp(-0.1*t);

% initial concentration
u_init = zeros(N,1);
u_init(1) = u_in_0;
u_init(2:end) = u_out;

%% main loop
% solution matrix
U = zeros(N,M);
U(:,1) = u_init;

for j = 2:M % time coordinate
    % system matrix
    Ku = zeros(N); % unconstrained system matrix
    for i = 2:N-1 % spatial (radial) coordinate
        Ku(i,i-1) = D - dr/2/r(i);
        Ku(i,i) = -(2*D + dr^2/dr);
        Ku(i,i+1) = D + dr/2/r(i);
    end
    K = Ku(2:N-1,2:N-1); % constrained matrix (after B.C. considered)
    
    % forcing vector
    fu = -dr^2/dt*U(:,j-1); % unconstrained forcing vector
    f = fu(2:N-1);
    
    % solution
    x_1 = u_in(j);
    x_end = u_out;
    x = K\(f - x_1*Ku(2:N-1,1) - x_end*Ku(2:N-1,N));
    
    % total solution
    U(:,j) = [x_1;x;x_end];
end
    
%% plots
figure;
idx_t = [1,2,3,5,10,50,M];
for ig = 1:length(idx_t)
    plot(r,U(:,idx_t(ig)),'LineWidth',2)
    hold on
end
hold off
xlabel('radial coordinate')
ylabel('concentration')


    
        
        
        
        
        
