function cost = cost_fn_kc_analysis_3lengths(B,Binv,kc_init1, w_init1, L1,  L2, L3, L4, p)
% function to generate the cost for the 3 length scenario
% for the kappa_c analysis
% the third length will be 90 mm, which will be read in
% 
% kc_init1 - The 90mm kappa_c 
% w_init1  - the initial rotation insertion
% L1       - the 90 mm insertion length
% L2, L3   - the two other insertion lengths of interest
% p        - the exponent of the length term
% 
% - written by: Dimitri Lezcano

%% Instantiations
    lengths = [L1 L2 L3 L4];

    kc_init2 = kappa_c_p(kc_init1,L1,L2,p);
    kc_init3 = kappa_c_p(kc_init1,L1,L3,p);
    kc_init4 = kappa_c_p(kc_init1,L1,L4,p);
%     r_90 = readmatrix("C:/Users/dlezcan1/Documents/Needle Shape Model/Kim Work/Dimitri/Position_Data_90mm.csv");
    
%     ds = L1 / (length(r_90)-1); % mm
    ds = 0.5;
    s = 0:ds:max([L1 L2 L3 L4]);
    
%% wv instantations for EP equation
    s1 = s(s <= L1);
    s2 = s(s <= L2);
    s3 = s(s <= L3);
    s4 = s(s <= L4);
    
    
    k01 = kc_init1*(1 - s1/L1).^2;
    k0_prime1 = -2*kc_init1/L1*(1 - s1/L1);
    w01 = [k01; zeros(1,length(s1)); zeros(1,length(s1))];
    w0_prime1 = [k0_prime1; zeros(1,length(s1)); zeros(1,length(s1))];
    
    k02 = kc_init2*(1 - s2/L2).^2;
    k0_prime2 = -2*kc_init2/L2*(1 - s2/L2);
    w02 = [k02; zeros(1,length(s2)); zeros(1,length(s2))];
    w0_prime2 = [k0_prime2; zeros(1,length(s2)); zeros(1,length(s2))];
    
    k03 = kc_init3*(1 - s3/L3).^2;
    k0_prime3 = -2*kc_init3/L3*(1 - s3/L3);
    w03 = [k03; zeros(1,length(s3)); zeros(1,length(s3))];
    w0_prime3 = [k0_prime3; zeros(1,length(s3)); zeros(1,length(s3))];
    
    k04 = kc_init4*(1 - s4/L4).^2;
    k0_prime4 = -2*kc_init4/L4*(1 - s4/L4);
    w04 = [k04; zeros(2,length(s4))];
    w0_prime4 = [k0_prime4; zeros(2,length(s4))];
    
    
%% EP calculations
    wv_1 = fn_intg_EP_Dimitri(w_init1, w01, w0_prime1, L1, s1, ds, B, Binv);
    wv_2 = fn_intg_EP_Dimitri(w_init1, w02, w0_prime2, L2, s2, ds, B, Binv);
    wv_3 = fn_intg_EP_Dimitri(w_init1, w03, w0_prime3, L3, s3, ds, B, Binv);
    wv_4 = fn_intg_EP_Dimitri(w_init1, w04, w0_prime4, L4, s4, ds, B, Binv);
    
    p1 = wv2r(wv_1, L1);
    p2 = wv2r(wv_2, L2);
    p3 = wv2r(wv_3, L3);
    p4 = wv2r(wv_4, L4);
    positions = {p1, p2, p3, p4};
    
%% Error Calculations
    cost = 0;
    for i = 1:4
        for j = i:4
            pi = positions{i};
            pj = positions{j};
            cost = cost + area_between_curves(pi, pj);
        end
    end
    
end


function cost = area_between_curves(pi, pj)
    % set-up for interpolation
    zi = pi(3,:);
    zj = pj(3,:);
    
    z_min = min(max(zi), max(zj));
    
    z = 0:.5:z_min;
    
    xi = interp1(zi, pi(1,:), z, 'spline');
    xj = interp1(zj, pj(1,:), z, 'spline');
    
    yi = interp1(zi, pi(2,:), z, 'spline');
    yj = interp1(zj, pj(2,:), z, 'spline');
    
    % take the area between curves.
    cost = sum(vecnorm([xi-xj; yi-yj]));
    
end
    
    
    
    