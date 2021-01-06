%% singlebend_needleshape
% this is a function to produce the singlebend needle shape
% given the measured curvatures of the needle
%
% - written by: Dimitri Lezcano

function pos = singlebend_needleshape(curvatures, aa_tip_locs, L)
% Input:
%   - curvatures: list of x-y curvatures measured at each of the AA locations
%           ( a #AA x 2 matrix )
%   - aa_locs: list of the active area locations (measured from the tip of
%   the needle)
%   - L: the needle length
%
    %% material properties
    % Stainless Steel 304
    Emod = 200e9*1e-6; % 200 GPa, conversion from N/m^2 to N/mm^2
    Pratio = 0.29; % Poisson's ratio
    diam = 0.9; % in mm
    Ibend = pi*diam^4/64;

    Gmod = Emod/2/(1+Pratio);
    Jtor = pi*diam^4/32;

    BendStiff = Emod*Ibend;
    TorStiff = Gmod*Jtor;

    B = diag([BendStiff,BendStiff,TorStiff]);
    Binv = inv(B);

    
    %% Needle arclength set-up
    aa_base_locs = L - aa_tip_locs;
    ds = 0.5;
    s = 0:ds:L;
    
    % get the arclength indices that are valid
    s_idx_aa = [];
    curvs_aa = [];
    for i = 1:length(aa_base_locs)
        aa_loc = aa_base_locs(i);
        idx = find(s == aa_loc);
        if ~isempty(idx)
            s_idx_aa = [s_idx_aa; idx];
            curvs_aa = [curvs_aa; curvatures(i,:), 0];
        end
    end
    
    % grab each of the valid arclengths and curvatures
    s_aa = s(s_idx_aa);
%     curvs_aa = curvatures(s_idx_aa, :);
    curvs_aa = [curvs_aa; zeros(size(curvs_aa, 1), 1)];
    
    
    %% Determine w_init and kc from measured curvatures
    
    