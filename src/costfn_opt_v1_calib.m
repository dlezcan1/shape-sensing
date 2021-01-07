function y = costfn_opt_v1_calib(eta,data,s_m,ds,N,B,Binv,scalef)

N_meas = size(data,2);
s_index_meas = s_m;

kc = eta;
w_init = [kc;0;0];
% w_init = eta(1:3);
% kc = eta(4);
% %alpha = eta(5);

L = (N-1)*ds; % in mm
s = [0:ds:L];

% linear
k0 = kc*ones(size(s));%*(1 - alpha*s);
w0 = [k0;zeros(1,N);zeros(1,N)];

%k0prime = -alpha*kc;
w0prime = zeros(3,N);%[k0prime*ones(1,N);zeros(1,N);zeros(1,N)];

% % exponential
% k0 = kc*exp(-alpha*s);
% w0 = [k0;zeros(1,N);zeros(1,N)];
% 
% k0prime = -alpha*k0;
% w0prime = [k0prime;zeros(1,N);zeros(1,N)];

wv = fn_intgEP_v1(w_init,w0,w0prime,0,ds,N,B,Binv);

yv = wv(1:3,s_index_meas) - data;
y = norm(yv,'fro')^2*scalef;

% y1 = 0;
% for i = 1:3
%     y1 = y1 + norm(yv(:,i),2)^2;
% end

