function y = costfn_inextrod_w0(x,m_ext,f_ext,B,Binv,ds,N,scalef)

e3 = [0;0;1];

%if length(x) == 3
    kc = x(1); % initial curvature
    p = x(2); % initial normal force
    f = x(3); % initial stretch/compression (and/or friction force)
% elseif length(x) == 2
%     kc = x(1);
%     p = x(2);
%     f = 0;
% end

[wv,fv] = fn_inextrod_w0(kc,p,f,m_ext,f_ext,B,Binv,ds,N);

% configuration
[Rmat,pmat] = cal_Rtraj(wv,ds,0);

% force and moment density vector viewed from R(0) (R(0) = I)
f_0_vec = zeros(3,N);
m_0_vec = zeros(3,N);
for i = 1:N
    f_0_vec(:,i) = Rmat(:,:,i)*f_ext(:,i);
    m_0_vec(:,i) = cross(pmat(:,i),Rmat(:,:,i)*f_ext(:,i)) + Rmat(:,:,i)*m_ext(:,i);
end
f_s = Simpson_vec_int(f_0_vec,ds);
m_s = Simpson_vec_int(m_0_vec,ds);

% cost function value
y1 = norm(fv(:,end),2) + 80*norm(wv(:,end),2);
y2 = norm(fv(:,1) - f_s,2) + norm(wv(:,1) - Binv*m_s,2);
y = y1 + y2;

y = y*scalef;


