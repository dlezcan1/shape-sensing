function y = costfn_opt_needle_doublebend_homlayer(eta,data,s_m,ix_s_turn,ds,N,B,Binv,scalef)

N_meas = size(data,2);
s_index_meas = s_m;

w_init = eta(1:3);
kc = eta(4);

L = (N-1)*ds; % in mm
s = [0:ds:L];

s1 = s(1:ix_s_turn);
s2 = s(ix_s_turn:end);
N1 = length(s1);
N2 = length(s2);

kc1 = kc*((s1(end) - s1(1))/L)^(2/3);
kc2 = kc*((s2(end) - s2(1))/L)^(2/3);

% intrinsic curvature kappa_0 (quadratic)
k0_1 = kc1*(1 - s1/L).^2;
% k0_1 = kc1*(s1(end) - s1).^2/L^2;
k0_2 = -kc2*(1 - s2/L).^2;
% k0_1 = kc1*(s1(end) - s1).^2/L^2;
% % k0_1 = kc1*(1 - s1/s1(end)).^2;
% k0_2 = -kc2*(s2 - s1(end)).^2/L^2.*(1 - s2/L);
% % k0_2 = -kc2*(1 - s2/s1(end)).^2.*(1 - s2/L);
k0_turn = (k0_1(end) + k0_2(1))/2; %0;
k0 = [k0_1(1:end-1),k0_turn,k0_2(2:end)];

k0prime1 = -2*kc1/L*(1 - s1/L);
% k0prime1 = -2*kc1/L^2*(s1(end) - s1);
k0prime2 = 2*kc2/L*(1 - s2/L);
% k0prime1 = -2*kc1/L^2*(s1(end) - s1);
% % k0prime1 = -2*kc1/s1(end)*(1 - s1);
% k0prime2 = -2*kc2/L^2*(s2 - s1(end)).*(1 - s2/L) + kc2/L^2*(s2 - s1(end)).^2/L; 
% % k0prime2 = -kc2/L.*(1 - s2/L) + kc2/L^2*(s2 - s1(end));
% % k0prime2 = 2*kc2/s1(end)*(1 - s2/s1(end)).*(1 - s2/L) + kc2/L*(1 - s2/s1(end)).^2;
k0prime_peak = (k0_2(2) - k0_1(end-1))/2/ds; 
k0prime = [k0prime1(1:end-1),k0prime_peak,k0prime2(2:end)];

% k0_1 = kc1*(s1(end)-s1).^2/L^2 - kc2*(1-s1(end)/L).*(1+s1(end)/L-2*s1/L);
% k0_2 = -kc2*(1-s2/L).^2;
% k0 = [k0_1(1:end-1),k0_2];
% 
% k0prime1 = -2*kc1/L^2*(s1(end) - s1) + kc2*(1-s1(end)/L).*(2/L);
% k0prime2 = 2*kc2/L*(1 - s2/L);
% k0prime = [k0prime1(1:end-1),k0prime2];

% intrinsic curvature \omega_0
w0 = [k0;zeros(size(s));zeros(size(s))];
w0prime = [k0prime;zeros(size(s));zeros(size(s))];
wv = fn_intgEP_needle(w_init,w0,w0prime,ds,N,B,Binv);
% w0_1 = [k0_1;zeros(size(s1));zeros(size(s1))];
% w0_2 = [k0_2;zeros(size(s2));zeros(size(s2))];
% 
% w0prime1 = [k0prime1;zeros(size(s1));zeros(size(s1))];
% w0prime2 = [k0prime2;zeros(size(s2));zeros(size(s2))];
% 
% wv_1 = fn_intgEP_v1(w_init,w0_1,w0prime1,0,ds,N1,B,Binv);
% wv_2 = fn_intgEP_v1(-wv_1(:,end),w0_2,w0prime2,s2(1),ds,N2,B,Binv);
% 
% wv = [wv_1(:,1:end-1),wv_2];

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

