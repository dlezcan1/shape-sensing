function y = costfn_opt_v2_singlebend(eta,data,s_m,ds,N,B,Binv,scalef)

N_meas = size(data,2);
s_index_meas = s_m;

w_init = eta(1:3);
kc = eta(4);

L = (N-1)*ds; % in mm
s = [0:ds:L];

% % intrinsic curvature (constant) --> not good
% k0 = kc*ones(1,N);
% w0 = [k0;zeros(2,N)];
%     
% w0prime = zeros(3,N);

% % intrinsic curvature (linear)
% k0 = kc*(1 - s/L);
% w0 = [k0;zeros(1,N);zeros(1,N)];
% 
% k0prime = -kc/L*ones(1,N);
% w0prime = [k0prime;zeros(1,N);zeros(1,N)];

% intrinsic curvature (quadratic)
k0 = kc*(1 - s/L).^2;
w0 = [k0;zeros(1,N);zeros(1,N)];

k0prime = -2*kc/L*(1 - s/L);
w0prime = [k0prime;zeros(1,N);zeros(1,N)];

% % intrinsic curvature (exponential)
% k0 = kc*exp(-s/L);
% w0 = [k0;zeros(1,N);zeros(1,N)];
% 
% k0prime = -kc/L*exp(-s/L);
% w0prime = [k0prime;zeros(1,N);zeros(1,N)];

% integration of the E-P equation
wv = fn_intgEP_v1(w_init,w0,w0prime,0,ds,N,B,Binv);

% % include zero torsion
% yv = wv(1:3,s_index_meas) - data;
% y = norm(yv,'fro')^2*scalef;

% exclude torsion
yv = wv(1:2,s_index_meas) - data(1:2,:);
y = norm(yv,'fro')^2*scalef;


% y1 = 0;
% for i = 1:3
%     y1 = y1 + norm(yv(:,i),2)^2;
% end

